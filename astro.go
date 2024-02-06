// Copyright 2024 Matthew P. Dargan. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Astro prints astronomical information.
//
// Usage:
//
//	astro [-jpokm] [-c nperiod] [-C tperiod] [-d date] [-e obj1 obj2] [-l nlat wlong elev] [-t ΔT]
//
// Astro reports upcoming celestial events, by default for 24 hours starting
// now.
//
// The -j flag causes astro to print the Julian date.
//
// The -p flag causes astro to print the positions of objects at the given time
// rather than searching for interesting conjunctions. For each, the name is
// followed by the right ascension (hours, minutes, seconds), declination
// (degrees, minutes, seconds), azimuth (degrees), elevation (degrees), and
// semidiameter (arc seconds). For the sun and moon, the magnitude is also
// printed. The first line of output presents the date and time, sidereal time,
// and the latitude, longitude, and elevation.
//
// The -o flag causes astro to search for stellar occultations.
//
// The -k flag causes astro to print times in local time (“kitchen clock”).
//
// The -m flag causes astro to include a single comet in the list of objects.
// This is modified (in the source) to refer to an approaching comet but in
// steady state usually refers to the last interesting comet (currently
// 153P/Ikeya–Zhang).
//
// The -c flag causes astro to report for n (default 1) successive days.
//
// The -C flag is used with -c and sets the interval to d days (or fractions of
// days).
//
// The -d flag causes astro to read the starting date.
//
// The -e flag causes astro to report distance between the centers of objects,
// in arc seconds, during eclipses or occultations involving obj1 and obj2.
//
// The -l flag causes astro to read the north latitude, west longitude, and
// elevation of the observation point. If l is missing, the initial position is
// read from the file $PLAN9/sky/here, or /usr/local/plan9/sky/here if $PLAN9
// is not set.
//
// The -t flag causes astro to read ΔT. ΔT is the difference between
// ephemeris and universal time (seconds) due to the slowing of the earth’s
// rotation. ΔT is normally calculated from an empirical formula. This option is
// needed only for very accurate timing of occultations, eclipses, etc.
package main

import (
	"bufio"
	"cmp"
	"errors"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"slices"
	"strconv"
	"strings"
	"time"
)

const (
	secondsPerDay       = 24 * 60 * 60
	meanSolarDaySeconds = 0.001704
	pipi                = 2 * math.Pi
	radian              = math.Pi / 180
	radsec              = radian / 3600
	converge            = 1e-14
	metersToFeet        = 3.28084
	per                 = 1.0
	npts                = 12
	deld                = per / npts
	maxe                = 0.999 // Can't do hyperbolas
	nevents             = 100
	dark                = 1 << iota
	signif
	ptime
	light
)

type obj1 struct {
	ra, decl2, semi2, az, el, mag float64
}

type obj2 struct {
	name, fname string
	f           func()
	point       [npts + 2]obj1
}

type obj3 struct {
	t1, e1, t2, e2, t3, e3, t4, e4, t5, e5 float64
}

type moonTab struct {
	f float64
	c [4]int
}

type evt struct {
	s    string
	tim  float64
	flag int
}

type occt struct {
	act, del0, del1, del2 obj1
}

var (
	julian    = flag.Bool("j", false, "print Julian date")
	pflag     = flag.Bool("p", false, "print positions of objects at the given time")
	oflag     = flag.Bool("o", false, "search for stellar occultations")
	local     = flag.Bool("k", false, "print times in local time")
	mflag     = flag.Bool("m", false, "include single comet in the list of objects")
	nperiods  = flag.Int("c", 1, "report for n successive days")
	p         = flag.Float64("C", per, "used with -c, set the interval to d days")
	startDate = flag.String("d", "", "read start date")
	eclipse   = flag.String("e", "", "report distance between the centers of objects")
	initLoc   = flag.String("l", "", "read latitude, longitude, and elevation")
	dt        = flag.Float64("t", 0.0, "read ΔT")

	root string
	wlong, awlong, nlat, elev,
	obliq, phi, eps, tobliq,
	dphi, deps,
	day, eday, capt, capt2, capt3, gst,
	ΔT, erad, glat,
	xms, yms, zms,
	xdot, ydot, zdot,
	ecc, incl, node, argp, mrad, anom, motion,
	lambda, beta, rad, mag, semi,
	alpha, delta, rp, hp,
	ra, decl, semi2,
	lha, decl2, lmb2,
	az, el,
	meday, seday, mhp, salph, sdelt, srad float64
	sao    string
	events []evt
	osun   = obj2{name: "sun", fname: "The sun", f: fsun}
	omoon  = obj2{name: "moon", fname: "The moon", f: moon}
	oshad  = obj2{name: "shadow", fname: "The shadow", f: shad}
	omerc  = obj2{name: "mercury", fname: "Mercury", f: merc}
	ovenus = obj2{name: "venus", fname: "Venus", f: venus}
	objs   = []*obj2{
		&osun, &omoon, &oshad, &omerc, &ovenus,
		{name: "mars", fname: "Mars", f: mars},
		{name: "jupiter", fname: "Jupiter", f: jup},
		{name: "saturn", fname: "Saturn", f: sat},
		{name: "uranus", fname: "Uranus", f: uran},
		{name: "neptune", fname: "Neptune", f: nept},
		{name: "pluto", fname: "Pluto", f: plut},
		{name: "comet", fname: "Comet", f: comet},
	}
	ostar      obj2
	occ        obj3
	occ1, occ2 occt
	moontab    = []moonTab{
		{f: 0.127, c: [4]int{0, 0, 0, 6}},
		{f: 13.902, c: [4]int{0, 0, 0, 4}},
		{f: 2369.912, c: [4]int{0, 0, 0, 2}},
		{f: 1.979, c: [4]int{1, 0, 0, 4}},
		{f: 191.953, c: [4]int{1, 0, 0, 2}},
		{f: 22639.500, c: [4]int{1, 0, 0, 0}},
		{f: -4586.465, c: [4]int{1, 0, 0, -2}},
		{f: -38.428, c: [4]int{1, 0, 0, -4}},
		{f: -.393, c: [4]int{1, 0, 0, -6}},
		{f: -.289, c: [4]int{0, 1, 0, 4}},
		{f: -24.420, c: [4]int{0, 1, 0, 2}},
		{f: -668.146, c: [4]int{0, 1, 0, 0}},
		{f: -165.145, c: [4]int{0, 1, 0, -2}},
		{f: -1.877, c: [4]int{0, 1, 0, -4}},
		{f: 0.403, c: [4]int{0, 0, 0, 3}},
		{f: -125.154, c: [4]int{0, 0, 0, 1}},
		{f: 0.213, c: [4]int{2, 0, 0, 4}},
		{f: 14.387, c: [4]int{2, 0, 0, 2}},
		{f: 769.016, c: [4]int{2, 0, 0, 0}},
		{f: -211.656, c: [4]int{2, 0, 0, -2}},
		{f: -30.773, c: [4]int{2, 0, 0, -4}},
		{f: -.570, c: [4]int{2, 0, 0, -6}},
		{f: -2.921, c: [4]int{1, 1, 0, 2}},
		{f: -109.673, c: [4]int{1, 1, 0, 0}},
		{f: -205.962, c: [4]int{1, 1, 0, -2}},
		{f: -4.391, c: [4]int{1, 1, 0, -4}},
		{f: -.072, c: [4]int{1, 1, 0, -6}},
		{f: 0.283, c: [4]int{1, -1, 0, 4}},
		{f: 14.577, c: [4]int{1, -1, 0, 2}},
		{f: 147.687, c: [4]int{1, -1, 0, 0}},
		{f: 28.475, c: [4]int{1, -1, 0, -2}},
		{f: 0.636, c: [4]int{1, -1, 0, -4}},
		{f: -.189, c: [4]int{0, 2, 0, 2}},
		{f: -7.486, c: [4]int{0, 2, 0, 0}},
		{f: -8.096, c: [4]int{0, 2, 0, -2}},
		{f: -.151, c: [4]int{0, 2, 0, -4}},
		{f: -.085, c: [4]int{0, 0, 2, 4}},
		{f: -5.741, c: [4]int{0, 0, 2, 2}},
		{f: -411.608, c: [4]int{0, 0, 2, 0}},
		{f: -55.173, c: [4]int{0, 0, 2, -2}},
		{f: -8.466, c: [4]int{1, 0, 0, 1}},
		{f: 18.609, c: [4]int{1, 0, 0, -1}},
		{f: 3.215, c: [4]int{1, 0, 0, -3}},
		{f: 0.150, c: [4]int{0, 1, 0, 3}},
		{f: 18.023, c: [4]int{0, 1, 0, 1}},
		{f: 0.560, c: [4]int{0, 1, 0, -1}},
		{f: 1.060, c: [4]int{3, 0, 0, 2}},
		{f: 36.124, c: [4]int{3, 0, 0, 0}},
		{f: -13.193, c: [4]int{3, 0, 0, -2}},
		{f: -1.187, c: [4]int{3, 0, 0, -4}},
		{f: -.293, c: [4]int{3, 0, 0, -6}},
		{f: -.290, c: [4]int{2, 1, 0, 2}},
		{f: -7.649, c: [4]int{2, 1, 0, 0}},
		{f: -8.627, c: [4]int{2, 1, 0, -2}},
		{f: -2.740, c: [4]int{2, 1, 0, -4}},
		{f: -.091, c: [4]int{2, 1, 0, -6}},
		{f: 1.181, c: [4]int{2, -1, 0, 2}},
		{f: 9.703, c: [4]int{2, -1, 0, 0}},
		{f: -2.494, c: [4]int{2, -1, 0, -2}},
		{f: 0.360, c: [4]int{2, -1, 0, -4}},
		{f: -1.167, c: [4]int{1, 2, 0, 0}},
		{f: -7.412, c: [4]int{1, 2, 0, -2}},
		{f: -.311, c: [4]int{1, 2, 0, -4}},
		{f: 0.757, c: [4]int{1, -2, 0, 2}},
		{f: 2.580, c: [4]int{1, -2, 0, 0}},
		{f: 2.533, c: [4]int{1, -2, 0, -2}},
		{f: -.103, c: [4]int{0, 3, 0, 0}},
		{f: -.344, c: [4]int{0, 3, 0, -2}},
		{f: -.992, c: [4]int{1, 0, 2, 2}},
		{f: -45.099, c: [4]int{1, 0, 2, 0}},
		{f: -.179, c: [4]int{1, 0, 2, -2}},
		{f: -.301, c: [4]int{1, 0, 2, -4}},
		{f: -6.382, c: [4]int{1, 0, -2, 2}},
		{f: 39.528, c: [4]int{1, 0, -2, 0}},
		{f: 9.366, c: [4]int{1, 0, -2, -2}},
		{f: 0.202, c: [4]int{1, 0, -2, -4}},
		{f: 0.415, c: [4]int{0, 1, 2, 0}},
		{f: -2.152, c: [4]int{0, 1, 2, -2}},
		{f: -1.440, c: [4]int{0, 1, -2, 2}},
		{f: 0.076, c: [4]int{0, 1, -2, 0}},
		{f: 0.384, c: [4]int{0, 1, -2, -2}},
		{f: -.586, c: [4]int{2, 0, 0, 1}},
		{f: 1.750, c: [4]int{2, 0, 0, -1}},
		{f: 1.225, c: [4]int{2, 0, 0, -3}},
		{f: 1.267, c: [4]int{1, 1, 0, 1}},
		{f: 0.137, c: [4]int{1, 1, 0, -1}},
		{f: 0.233, c: [4]int{1, 1, 0, -3}},
		{f: -.122, c: [4]int{1, -1, 0, 1}},
		{f: -1.089, c: [4]int{1, -1, 0, -1}},
		{f: -.276, c: [4]int{1, -1, 0, -3}},
		{f: 0.255, c: [4]int{0, 0, 2, 1}},
		{f: 0.584, c: [4]int{0, 0, 2, -1}},
		{f: 0.254, c: [4]int{0, 0, 2, -3}},
		{f: 0.070, c: [4]int{4, 0, 0, 2}},
		{f: 1.938, c: [4]int{4, 0, 0, 0}},
		{f: -.952, c: [4]int{4, 0, 0, -2}},
		{f: -.551, c: [4]int{3, 1, 0, 0}},
		{f: -.482, c: [4]int{3, 1, 0, -2}},
		{f: -.100, c: [4]int{3, 1, 0, -4}},
		{f: 0.088, c: [4]int{3, -1, 0, 2}},
		{f: 0.681, c: [4]int{3, -1, 0, 0}},
		{f: -.183, c: [4]int{3, -1, 0, -2}},
		{f: -.297, c: [4]int{2, 2, 0, -2}},
		{f: -.161, c: [4]int{2, 2, 0, -4}},
		{f: 0.197, c: [4]int{2, -2, 0, 0}},
		{f: 0.254, c: [4]int{2, -2, 0, -2}},
		{f: -.250, c: [4]int{1, 3, 0, -2}},
		{f: -.123, c: [4]int{2, 0, 2, 2}},
		{f: -3.996, c: [4]int{2, 0, 2, 0}},
		{f: 0.557, c: [4]int{2, 0, 2, -2}},
		{f: -.459, c: [4]int{2, 0, -2, 2}},
		{f: -1.370, c: [4]int{2, 0, -2, 0}},
		{f: 0.538, c: [4]int{2, 0, -2, -2}},
		{f: 0.173, c: [4]int{2, 0, -2, -4}},
		{f: 0.263, c: [4]int{1, 1, 2, 0}},
		{f: 0.083, c: [4]int{1, 1, -2, 2}},
		{f: -.083, c: [4]int{1, 1, -2, 0}},
		{f: 0.426, c: [4]int{1, 1, -2, -2}},
		{f: -.304, c: [4]int{1, -1, 2, 0}},
		{f: -.372, c: [4]int{1, -1, -2, 2}},
		{f: 0.083, c: [4]int{1, -1, -2, 0}},
		{f: 0.418, c: [4]int{0, 0, 4, 0}},
		{f: 0.074, c: [4]int{0, 0, 4, -2}},
		{f: 0.130, c: [4]int{3, 0, 0, -1}},
		{f: 0.092, c: [4]int{2, 1, 0, 1}},
		{f: 0.084, c: [4]int{2, 1, 0, -3}},
		{f: -.352, c: [4]int{2, -1, 0, -1}},
		{f: 0.113, c: [4]int{5, 0, 0, 0}},
		{f: -.330, c: [4]int{3, 0, 2, 0}},
		{f: 0.090, c: [4]int{1, 0, 4, 0}},
		{f: -.080, c: [4]int{1, 0, -4, 0}},
		{},
		{f: -112.79, c: [4]int{0, 0, 0, 1}},
		{f: 2373.36, c: [4]int{0, 0, 0, 2}},
		{f: -4.01, c: [4]int{0, 0, 0, 3}},
		{f: 14.06, c: [4]int{0, 0, 0, 4}},
		{f: 6.98, c: [4]int{1, 0, 0, 4}},
		{f: 192.72, c: [4]int{1, 0, 0, 2}},
		{f: -13.51, c: [4]int{1, 0, 0, 1}},
		{f: 22609.07, c: [4]int{1, 0, 0, 0}},
		{f: 3.59, c: [4]int{1, 0, 0, -1}},
		{f: -4578.13, c: [4]int{1, 0, 0, -2}},
		{f: 5.44, c: [4]int{1, 0, 0, -3}},
		{f: -38.64, c: [4]int{1, 0, 0, -4}},
		{f: 14.78, c: [4]int{2, 0, 0, 2}},
		{f: 767.96, c: [4]int{2, 0, 0, 0}},
		{f: 2.01, c: [4]int{2, 0, 0, -1}},
		{f: -152.53, c: [4]int{2, 0, 0, -2}},
		{f: -34.07, c: [4]int{2, 0, 0, -4}},
		{f: 2.96, c: [4]int{3, 0, 0, 2}},
		{f: 50.64, c: [4]int{3, 0, 0, 0}},
		{f: -16.40, c: [4]int{3, 0, 0, -2}},
		{f: 3.60, c: [4]int{4, 0, 0, 0}},
		{f: -1.58, c: [4]int{4, 0, 0, -2}},
		{f: -1.59, c: [4]int{0, 1, 0, 4}},
		{f: -25.10, c: [4]int{0, 1, 0, 2}},
		{f: 17.93, c: [4]int{0, 1, 0, 1}},
		{f: -126.98, c: [4]int{0, 1, 0, 0}},
		{f: -165.06, c: [4]int{0, 1, 0, -2}},
		{f: -6.46, c: [4]int{0, 1, 0, -4}},
		{f: -1.68, c: [4]int{0, 2, 0, 2}},
		{f: -16.35, c: [4]int{0, 2, 0, -2}},
		{f: -11.75, c: [4]int{1, 1, 0, 2}},
		{f: 1.52, c: [4]int{1, 1, 0, 1}},
		{f: -115.18, c: [4]int{1, 1, 0, 0}},
		{f: -182.36, c: [4]int{1, 1, 0, -2}},
		{f: -9.66, c: [4]int{1, 1, 0, -4}},
		{f: -2.27, c: [4]int{-1, 1, 0, 4}},
		{f: -23.59, c: [4]int{-1, 1, 0, 2}},
		{f: -138.76, c: [4]int{-1, 1, 0, 0}},
		{f: -31.70, c: [4]int{-1, 1, 0, -2}},
		{f: -1.53, c: [4]int{-1, 1, 0, -4}},
		{f: -10.56, c: [4]int{2, 1, 0, 0}},
		{f: -7.59, c: [4]int{2, 1, 0, -2}},
		{f: -2.54, c: [4]int{2, 1, 0, -4}},
		{f: 3.32, c: [4]int{2, -1, 0, 2}},
		{f: 11.67, c: [4]int{2, -1, 0, 0}},
		{f: -6.12, c: [4]int{1, 2, 0, -2}},
		{f: -2.40, c: [4]int{-1, 2, 0, 2}},
		{f: -2.32, c: [4]int{-1, 2, 0, 0}},
		{f: -1.82, c: [4]int{-1, 2, 0, -2}},
		{f: -52.14, c: [4]int{0, 0, 2, -2}},
		{f: -1.67, c: [4]int{0, 0, 2, -4}},
		{f: -9.52, c: [4]int{1, 0, 2, -2}},
		{f: -85.13, c: [4]int{-1, 0, 2, 0}},
		{f: 3.37, c: [4]int{-1, 0, 2, -2}},
		{f: -2.26, c: [4]int{0, 1, 2, -2}},
		{},
		{f: -0.725, c: [4]int{0, 0, 0, 1}},
		{f: 0.601, c: [4]int{0, 0, 0, 2}},
		{f: 0.394, c: [4]int{0, 0, 0, 3}},
		{f: -.445, c: [4]int{1, 0, 0, 4}},
		{f: 0.455, c: [4]int{1, 0, 0, 1}},
		{f: 0.192, c: [4]int{1, 0, 0, -3}},
		{f: 5.679, c: [4]int{2, 0, 0, -2}},
		{f: -.308, c: [4]int{2, 0, 0, -4}},
		{f: -.166, c: [4]int{3, 0, 0, 2}},
		{f: -1.300, c: [4]int{3, 0, 0, 0}},
		{f: 0.258, c: [4]int{3, 0, 0, -2}},
		{f: -1.302, c: [4]int{0, 1, 0, 0}},
		{f: -.416, c: [4]int{0, 1, 0, -4}},
		{f: -.740, c: [4]int{0, 2, 0, -2}},
		{f: 0.787, c: [4]int{1, 1, 0, 2}},
		{f: 0.461, c: [4]int{1, 1, 0, 0}},
		{f: 2.056, c: [4]int{1, 1, 0, -2}},
		{f: -.471, c: [4]int{1, 1, 0, -4}},
		{f: -.443, c: [4]int{-1, 1, 0, 2}},
		{f: 0.679, c: [4]int{-1, 1, 0, 0}},
		{f: -1.540, c: [4]int{-1, 1, 0, -2}},
		{f: 0.259, c: [4]int{2, 1, 0, 0}},
		{f: -.212, c: [4]int{2, -1, 0, 2}},
		{f: -.151, c: [4]int{2, -1, 0, 0}},
		{},
		{f: -526.069, c: [4]int{0, 0, 1, -2}},
		{f: -3.352, c: [4]int{0, 0, 1, -4}},
		{f: 44.297, c: [4]int{1, 0, 1, -2}},
		{f: -6.000, c: [4]int{1, 0, 1, -4}},
		{f: 20.599, c: [4]int{-1, 0, 1, 0}},
		{f: -30.598, c: [4]int{-1, 0, 1, -2}},
		{f: -24.649, c: [4]int{-2, 0, 1, 0}},
		{f: -2.000, c: [4]int{-2, 0, 1, -2}},
		{f: -22.571, c: [4]int{0, 1, 1, -2}},
		{f: 10.985, c: [4]int{0, -1, 1, -2}},
		{},
		{f: 0.2607, c: [4]int{0, 0, 0, 4}},
		{f: 28.2333, c: [4]int{0, 0, 0, 2}},
		{f: 0.0433, c: [4]int{1, 0, 0, 4}},
		{f: 3.0861, c: [4]int{1, 0, 0, 2}},
		{f: 186.5398, c: [4]int{1, 0, 0, 0}},
		{f: 34.3117, c: [4]int{1, 0, 0, -2}},
		{f: 0.6008, c: [4]int{1, 0, 0, -4}},
		{f: -.3000, c: [4]int{0, 1, 0, 2}},
		{f: -.3997, c: [4]int{0, 1, 0, 0}},
		{f: 1.9178, c: [4]int{0, 1, 0, -2}},
		{f: 0.0339, c: [4]int{0, 1, 0, -4}},
		{f: -0.9781, c: [4]int{0, 0, 0, 1}},
		{f: 0.2833, c: [4]int{2, 0, 0, 2}},
		{f: 10.1657, c: [4]int{2, 0, 0, 0}},
		{f: -.3039, c: [4]int{2, 0, 0, -2}},
		{f: 0.3722, c: [4]int{2, 0, 0, -4}},
		{f: 0.0109, c: [4]int{2, 0, 0, -6}},
		{f: -.0484, c: [4]int{1, 1, 0, 2}},
		{f: -.9490, c: [4]int{1, 1, 0, 0}},
		{f: 1.4437, c: [4]int{1, 1, 0, -2}},
		{f: 0.0673, c: [4]int{1, 1, 0, -4}},
		{f: 0.2302, c: [4]int{1, -1, 0, 2}},
		{f: 1.1528, c: [4]int{1, -1, 0, 0}},
		{f: -.2257, c: [4]int{1, -1, 0, -2}},
		{f: -.0102, c: [4]int{1, -1, 0, -4}},
		{f: 0.0918, c: [4]int{0, 2, 0, -2}},
		{f: -.0124, c: [4]int{0, 0, 2, 0}},
		{f: -.1052, c: [4]int{0, 0, 2, -2}},
		{f: -.1093, c: [4]int{1, 0, 0, 1}},
		{f: 0.0118, c: [4]int{1, 0, 0, -1}},
		{f: -.0386, c: [4]int{1, 0, 0, -3}},
		{f: 0.1494, c: [4]int{0, 1, 0, 1}},
		{f: 0.0243, c: [4]int{3, 0, 0, 2}},
		{f: 0.6215, c: [4]int{3, 0, 0, 0}},
		{f: -.1187, c: [4]int{3, 0, 0, -2}},
		{f: -.1038, c: [4]int{2, 1, 0, 0}},
		{f: -.0192, c: [4]int{2, 1, 0, -2}},
		{f: 0.0324, c: [4]int{2, 1, 0, -4}},
		{f: 0.0213, c: [4]int{2, -1, 0, 2}},
		{f: 0.1268, c: [4]int{2, -1, 0, 0}},
		{f: -.0106, c: [4]int{1, 2, 0, 0}},
		{f: 0.0484, c: [4]int{1, 2, 0, -2}},
		{f: 0.0112, c: [4]int{1, -2, 0, 2}},
		{f: 0.0196, c: [4]int{1, -2, 0, 0}},
		{f: -.0212, c: [4]int{1, -2, 0, -2}},
		{f: -.0833, c: [4]int{1, 0, 2, -2}},
		{f: -.0481, c: [4]int{1, 0, -2, 2}},
		{f: -.7136, c: [4]int{1, 0, -2, 0}},
		{f: -.0112, c: [4]int{1, 0, -2, -2}},
		{f: -.0100, c: [4]int{2, 0, 0, 1}},
		{f: 0.0155, c: [4]int{2, 0, 0, -1}},
		{f: 0.0164, c: [4]int{1, 1, 0, 1}},
		{f: 0.0401, c: [4]int{4, 0, 0, 0}},
		{f: -.0130, c: [4]int{4, 0, 0, -2}},
		{f: 0.0115, c: [4]int{3, -1, 0, 0}},
		{f: -.0141, c: [4]int{2, 0, -2, -2}},
		{},
	}
	sunf = [8][]float64{
		{
			-.265, 0,
			3.760, 0,
			.200, 0,
		},
		{
			-.021, 0,
			5.180, 0,
			1.882, 3.8991,
			-.030, 0,
		},
		{
			0.075, 5.1766,
			4.838, 5.2203,
			0.074, 3.6285,
			0.116, 2.5988,
			5.526, 2.5885,
			2.497, 5.5143,
			0.044, 5.4350,
			0.666, 3.1016,
			1.559, 6.0258,
			1.024, 5.5527,
			0.210, 3.5989,
			0.144, 3.4104,
			0.152, 6.0004,
			0.084, 4.1120,
			0.037, 3.8711,
			0.123, 3.4086,
			0.154, 6.2762,
			0.038, 4.6094,
			0.020, 5.1313,
			0.042, 4.5239,
			0.032, 0.8517,
			0.273, 3.7996,
			0.048, 4.5431,
			0.041, 6.0388,
			2.043, 6.0020,
			1.770, 3.4977,
			0.028, 2.5831,
			0.129, 5.1348,
			0.425, 5.9146,
			0.034, 1.2391,
			0.500, 1.8357,
			0.585, 5.8304,
			0.085, 0.9529,
			0.204, 1.7593,
			0.020, 3.2463,
			0.154, 3.9689,
			0.101, 1.6808,
			0.049, 3.0805,
			0.106, 3.8868,
			0.052, 6.0895,
			0.021, 3.7559,
			0.028, 5.2011,
			0.062, 6.0388,
			0.044, 1.8483,
			0.045, 3.9759,
			0.021, 5.3931,
			0.026, 1.9722,
			0.163, 3.4662,
			7.208, 3.1334,
			2.600, 4.5940,
			0.073, 4.8223,
			0.069, 1.4102,
			2.731, 1.5210,
			1.610, 1.9110,
			0.073, 4.4087,
			0.164, 2.9758,
			0.556, 1.4425,
			0.210, 1.7261,
			0.044, 2.9356,
			0.080, 1.3561,
			0.419, 1.7555,
			0.320, 4.7030,
			0.108, 5.0719,
			0.112, 5.1243,
			0.021, 5.0440,
		},
		{
			6.454, 0,
			0.177, 0,
			-.424, 0,
			0.039, 0,
			-.064, 0,
			0.172, 0,
		},
		{
			-.092, 1.6354,
			-.067, 2.1468,
			-.210, 2.6494,
			-.166, 4.6338,
		},
		{
			0.576, 0,
			-.047, 0,
			0.021, 0,
		},
		{
			2.359e-6, 3.6607,
			6.842e-6, 1.0180,
			0.869e-6, 3.9567,
			1.045e-6, 1.5332,
			1.497e-6, 4.4691,
			0.376e-6, 2.0295,
			2.057e-6, 4.0941,
			0.215e-6, 4.3459,
			0.478e-6, 0.2648,
			0.208e-6, 1.9548,
			7.067e-6, 1.5630,
			0.244e-6, 5.9097,
			4.026e-6, 6.2526,
			1.459e-6, 0.3409,
			0.281e-6, 1.4172,
			0.803e-6, 6.1533,
			0.429e-6, 0.1850,
		},
		{
			13.36e-6, 0,
			-1.33e-6, 0,
			0.37e-6, 0,
			0.36e-6, 0,
		},
	}
	sunc = [8][]int{
		{
			4, -7, 3, 0,
			-8, 4, 0, 3,
			15, -8, 0, 0,
		},
		{
			4, -7, 3, 0, 0,
			-8, 4, 0, 3, 0,
			0, 13, -8, 0, 1,
			15, -8, 0, 0, 0,
		},
		{
			0, 0, 1, 0, 0,
			0, -1, 1, 0, 0,
			0, -2, 1, 0, 0,
			0, -1, 2, 0, 0,
			0, -2, 2, 0, 0,
			0, -3, 2, 0, 0,
			0, -4, 2, 0, 0,
			0, -3, 3, 0, 0,
			0, -4, 3, 0, 0,
			0, -5, 3, 0, 0,
			0, -4, 4, 0, 0,
			0, -5, 4, 0, 0,
			0, -6, 4, 0, 0,
			0, -5, 5, 0, 0,
			0, -6, 5, 0, 0,
			0, -7, 5, 0, 0,
			0, -8, 5, 0, 0,
			0, -6, 6, 0, 0,
			0, -7, 7, 0, 0,
			0, -12, 8, 0, 0,
			0, -14, 8, 0, 0,
			-1, 1, 0, 0, 0,
			-1, 0, 0, 0, 0,
			-2, 3, 0, 0, 0,
			-2, 2, 0, 0, 0,
			-2, 1, 0, 0, 0,
			-2, 0, 0, 0, 0,
			-3, 3, 0, 0, 0,
			-3, 2, 0, 0, 0,
			-4, 4, 0, 0, 0,
			-4, 3, 0, 0, 0,
			-4, 2, 0, 0, 0,
			-5, 4, 0, 0, 0,
			-5, 3, 0, 0, 0,
			-6, 5, 0, 0, 0,
			-6, 4, 0, 0, 0,
			-6, 3, 0, 0, 0,
			-7, 5, 0, 0, 0,
			-7, 4, 0, 0, 0,
			-8, 5, 0, 0, 0,
			-8, 4, 0, 0, 0,
			-9, 6, 0, 0, 0,
			-9, 5, 0, 0, 0,
			-11, 6, 0, 0, 0,
			-13, 7, 0, 0, 0,
			-15, 9, 0, 0, 0,
			-17, 9, 0, 0, 0,
			0, 2, 0, -1, 0,
			0, 1, 0, -1, 0,
			0, 0, 0, -1, 0,
			0, -1, 0, -1, 0,
			0, 3, 0, -2, 0,
			0, 2, 0, -2, 0,
			0, 1, 0, -2, 0,
			0, 0, 0, -2, 0,
			0, 3, 0, -3, 0,
			0, 2, 0, -3, 0,
			0, 1, 0, -3, 0,
			0, 3, 0, -4, 0,
			0, 2, 0, -4, 0,
			0, 1, 0, 0, -1,
			0, 0, 0, 0, -1,
			0, 2, 0, 0, -2,
			0, 1, 0, 0, -2,
			0, 2, 0, 0, -3,
		},
		{
			1, 0, 0,
			1, 1, 0,
			1, -1, 0,
			3, -1, 0,
			1, 0, 1,
			1, 0, -1,
		},
		{
			-2, 1, 0,
			-3, 2, 0,
			-4, 2, 0,
			1, 0, -2,
		},
		{
			1, 0, 0,
			1, -1, 0,
			-1, 0, 2,
		},
		{
			0, -1, 1, 0, 0,
			0, -2, 2, 0, 0,
			0, -3, 2, 0, 0,
			0, -3, 3, 0, 0,
			0, -4, 3, 0, 0,
			0, -4, 4, 0, 0,
			-2, 2, 0, 0, 0,
			-3, 2, 0, 0, 0,
			-4, 3, 0, 0, 0,
			0, 2, 0, -1, 0,
			0, 1, 0, -1, 0,
			0, 0, 0, -1, 0,
			0, 2, 0, -2, 0,
			0, 1, 0, -2, 0,
			0, 3, 0, -3, 0,
			0, 2, 0, -3, 0,
			0, 1, 0, 0, -1,
		},
		{
			1, 0, 0,
			1, -1, 0,
			1, 1, 0,
			1, 0, -1,
		},
	}
	mercf = [7][]float64{
		{
			0.013, 0.6807,
			0.048, 0.6283,
			0.185, 0.6231,
			0.711, 0.6191,
			0.285, 0.5784,
			0.075, 0.5411,
			0.019, 0.5585,
			0.010, 2.8449,
			0.039, 2.8117,
			0.147, 2.8135,
			0.552, 2.8126,
			2.100, 2.8126,
			3.724, 2.8046,
			0.729, 2.7883,
			0.186, 2.7890,
			0.049, 2.7943,
			0.013, 2.7402,
			0.033, 1.8361,
			0.118, 1.8396,
			0.431, 1.8391,
			1.329, 1.8288,
			0.539, 4.8686,
			0.111, 4.8904,
			0.027, 4.8956,
			0.012, 3.9794,
			0.056, 3.9636,
			0.294, 3.9910,
			0.484, 3.9514,
			0.070, 3.9270,
			0.018, 3.9270,
			0.013, 6.1261,
			0.050, 6.1052,
			0.185, 6.1069,
			0.685, 6.1011,
			2.810, 6.1062,
			7.356, 6.0699,
			1.471, 6.0685,
			0.375, 6.0687,
			0.098, 6.0720,
			0.026, 6.0476,
			0.062, 5.1540,
			0.122, 5.1191,
			0.011, 0.9076,
			0.074, 1.0123,
			0.106, 0.9372,
			0.017, 0.9425,
			0.020, 0.0506,
			0.052, 0.0384,
			0.052, 3.0281,
			0.012, 3.0543,
			0.011, 2.1642,
			0.016, 2.2340,
			0.040, 4.3912,
			0.080, 4.4262,
			0.016, 4.4506,
		},
		{
			0.014, 1.0996,
			0.056, 1.1153,
			0.219, 1.1160,
			0.083, 1.0734,
			0.024, 0.9442,
			0.018, 3.8432,
			0.070, 3.8293,
			0.256, 3.8230,
			0.443, 3.8132,
			0.080, 3.7647,
			0.020, 3.7734,
			0.019, 0.0000,
			0.133, 0.1134,
			0.129, 6.2588,
			0.026, 6.2413,
			0.026, 2.6599,
			0.087, 2.6232,
			0.374, 2.6496,
			0.808, 2.5470,
			0.129, 2.5587,
			0.019, 2.5534,
			0.012, 2.1642,
		},
		{
			0.014, 3.1416,
			0.047, 3.1625,
			0.179, 3.1695,
			0.697, 3.1603,
			0.574, 4.1315,
			0.181, 4.2537,
			0.047, 4.2481,
			0.013, 4.2062,
			0.018, 0.6650,
			0.069, 0.6405,
			0.253, 0.6449,
			0.938, 0.6454,
			3.275, 0.6458,
			0.499, 0.5569,
			0.119, 0.5271,
			0.032, 0.5184,
			0.030, 0.4939,
			0.106, 0.4171,
			0.353, 0.4510,
			0.056, 0.3840,
			0.013, 0.3142,
			0.028, 0.2531,
		},
		{
			0.034, 0.9512,
			0.060, 4.7962,
			0.028, 4.7124,
			0.028, 4.1836,
			0.102, 4.1871,
			0.380, 4.1864,
			0.059, 4.1818,
			0.015, 4.2185,
			0.012, 4.1713,
			0.050, 4.1870,
		},
		{
			0.218e-6, 5.3369,
			0.491e-6, 5.3281,
			0.172e-6, 2.1642,
			0.091e-6, 2.1084,
			0.204e-6, 1.2460,
			0.712e-6, 1.2413,
			2.370e-6, 1.2425,
			0.899e-6, 1.2303,
			0.763e-6, 4.3633,
			0.236e-6, 4.3590,
			0.163e-6, 0.2705,
			0.541e-6, 0.2710,
			1.157e-6, 0.2590,
			0.099e-6, 0.1798,
			0.360e-6, 2.4237,
			0.234e-6, 2.3740,
			0.253e-6, 4.5365,
			0.849e-6, 4.5293,
			2.954e-6, 4.5364,
			0.282e-6, 4.4581,
			1.550e-6, 1.3570,
			0.472e-6, 1.3561,
			0.135e-6, 1.3579,
			0.081e-6, 3.5936,
			0.087e-6, 3.5500,
			0.087e-6, 5.7334,
		},
		{
			0.181e-6, 5.8275,
			0.095e-6, 2.2427,
			0.319e-6, 2.2534,
			0.256e-6, 2.2403,
			0.157e-6, 4.8292,
			0.106e-6, 1.0332,
			0.397e-6, 1.0756,
			0.143e-6, 4.0980,
		},
		{
			0.222e-6, 1.6024,
			0.708e-6, 1.5949,
			0.191e-6, 5.7914,
			0.100e-6, 5.3564,
			0.347e-6, 5.3548,
			1.185e-6, 5.3576,
			3.268e-6, 5.3579,
			0.371e-6, 2.2148,
			0.160e-6, 2.1241,
			0.134e-6, 5.1260,
			0.347e-6, 5.1620,
		},
	}
	mercc = [7][]int{
		{
			4, 1,
			3, 1,
			2, 1,
			1, 1,
			0, 1,
			-1, 1,
			-2, 1,
			6, 2,
			5, 2,
			4, 2,
			3, 2,
			2, 2,
			1, 2,
			0, 2,
			-1, 2,
			-2, 2,
			-3, 2,
			5, 3,
			4, 3,
			3, 3,
			2, 3,
			1, 3,
			0, 3,
			-1, 3,
			5, 4,
			4, 4,
			3, 4,
			2, 4,
			1, 4,
			0, 4,
			7, 5,
			6, 5,
			5, 5,
			4, 5,
			3, 5,
			2, 5,
			1, 5,
			0, 5,
			-1, 5,
			-2, 5,
			4, 6,
			3, 6,
			5, 7,
			4, 7,
			3, 7,
			2, 7,
			5, 8,
			4, 8,
			3, 8,
			2, 8,
			5, 9,
			4, 9,
			5, 10,
			4, 10,
			3, 10,
		},
		{
			3, 1,
			2, 1,
			1, 1,
			0, 1,
			-1, 1,
			4, 2,
			3, 2,
			2, 2,
			1, 2,
			0, 2,
			-1, 2,
			3, 3,
			2, 3,
			1, 3,
			0, 3,
			4, 4,
			3, 4,
			2, 4,
			1, 4,
			0, 4,
			-1, 4,
			2, 5,
		},
		{
			4, 1,
			3, 1,
			2, 1,
			1, 1,
			0, 1,
			-1, 1,
			-2, 1,
			-3, 1,
			5, 2,
			4, 2,
			3, 2,
			2, 2,
			1, 2,
			0, 2,
			-1, 2,
			-2, 2,
			3, 3,
			2, 3,
			1, 3,
			0, 3,
			-1, 3,
			1, 4,
		},
		{
			1, 1,
			0, 1,
			-1, 1,
			3, 2,
			2, 2,
			1, 2,
			0, 2,
			-1, 2,
			2, 3,
			1, 3,
		},
		{
			2, 1,
			1, 1,
			0, 1,
			-1, 1,
			4, 2,
			3, 2,
			2, 2,
			1, 2,
			0, 2,
			-1, 2,
			4, 3,
			3, 3,
			2, 3,
			0, 3,
			3, 4,
			2, 4,
			5, 5,
			4, 5,
			3, 5,
			2, 5,
			1, 5,
			0, 5,
			-1, 5,
			4, 6,
			3, 6,
			4, 7,
		},
		{
			1, 1,
			3, 2,
			2, 2,
			1, 2,
			2, 3,
			3, 4,
			2, 4,
			0, 4,
		},
		{
			2, 1,
			1, 1,
			-1, 1,
			4, 2,
			3, 2,
			2, 2,
			1, 2,
			0, 2,
			-1, 2,
			2, 3,
			1, 3,
		},
	}
	venf = [3][]float64{
		{
			4.889, 2.0788,
			11.261, 2.5870,
			7.128, 6.2384,
			3.446, 2.3721,
			1.034, 0.4632,
			1.575, 3.3847,
			1.439, 2.4099,
			1.208, 4.1464,
			2.966, 3.6318,
			1.563, 4.6829,
		},
		{
			0.122, 4.2726,
			0.300, 0.0218,
			0.159, 1.3491,
		},
		{
			2.246e-6, 0.5080,
			9.772e-6, 1.0159,
			8.271e-6, 4.6674,
			0.737e-6, 0.8267,
			1.426e-6, 5.1747,
			0.510e-6, 5.7009,
			1.572e-6, 1.8188,
			0.717e-6, 2.2969,
			2.991e-6, 2.0611,
			1.335e-6, 0.9628,
		},
	}
	venc = [3][]int{
		{
			1, -1, 0, 0,
			2, -2, 0, 0,
			3, -3, 0, 0,
			2, -3, 0, 0,
			4, -4, 0, 0,
			4, -5, 0, 0,
			3, -5, 0, 0,
			1, 0, -3, 0,
			1, 0, 0, -1,
			0, 0, 0, -1,
		},
		{
			0, -1, 0,
			4, -5, 0,
			1, 0, -2,
		},
		{
			1, -1, 0, 0,
			2, -2, 0, 0,
			3, -3, 0, 0,
			2, -3, 0, 0,
			4, -4, 0, 0,
			5, -5, 0, 0,
			4, -5, 0, 0,
			2, 0, -3, 0,
			1, 0, 0, -1,
			2, 0, 0, -2,
		},
	}
	nutf = [4][]float64{
		{
			.2088, 0,
			-1.2730, 0,
			.1258, 0,
			-.0496, 0,
			.0214, 0,
			.0124, 0,
		},
		{
			9.2109, 0,
			-.0904, 0,
			.5519, 0,
			.0215, 0,
			-.0093, 0,
			-.0066, 0,
		},
		{
			-.2037, 0,
			.0675, 0,
			-.0342, 0,
			-.0261, 0,
			-.0149, 0,
			.0114, 0,
			.0060, 0,
			.0058, 0,
			-.0057, 0,
			-.0052, 0,
		},
		{
			.0884, 0,
			.0183, 0,
			.0113, 0,
			-.0050, 0,
		},
	}
	nutc = [4][]int{
		{
			2, 0, 0, 0,
			2, 2, -2, 0,
			0, 0, 0, 1,
			2, 2, -2, 1,
			2, 2, -2, -1,
			1, 2, -2, 0,
		},
		{
			1, 0, 0, 0,
			2, 0, 0, 0,
			2, 2, -2, 0,
			2, 2, -2, 1,
			2, 2, -2, -1,
			1, 2, -2, 0,
		},
		{
			2, 2, 0, 0,
			0, 0, 1, 0,
			1, 2, 0, 0,
			2, 2, 1, 0,
			0, 0, 1, -2,
			2, 2, -1, 0,
			0, 0, 0, 2,
			1, 0, 1, 0,
			1, 0, -1, 0,
			2, 2, -1, 2,
		},
		{
			2, 2, 0,
			1, 2, 0,
			2, 2, 1,
			2, 2, -1,
		},
	}
	elemUran = []float64{
		36525.0,     // eday of epoc
		19.19126393, // semi major axis (au)
		0.04716771,  // eccentricity
		0.76986,     // inclination (deg)
		74.22988,    // longitude of the ascending node (deg)
		170.96424,   // longitude of perihelion (deg)
		313.23218,   // mean longitude (deg)
		0.00152025,  // (au/julian century)
		-0.00019150, // (e/julian century)
		-2.09,       // (arcsec/julian century)
		-1681.40,    // (arcsec/julian century)
		1312.56,     // (arcsec/julian century)
		1542547.79,  // (arcsec/julian century)
	}
	elemNept = []float64{
		36525.0,     // eday of epoc
		30.06896348, // semi major axis (au)
		0.00858587,  // eccentricity
		1.76917,     // inclination (deg)
		131.72169,   // longitude of the ascending node (deg)
		44.97135,    // longitude of perihelion (deg)
		304.88003,   // mean longitude (deg)
		-0.00125196, // (au/julian century)
		0.0000251,   // (e/julian century)
		-3.64,       // (arcsec/julian century)
		-151.25,     // (arcsec/julian century)
		-844.43,     // (arcsec/julian century)
		786449.21,   // (arcsec/julian century)
	}
	elemPlut = []float64{
		36525.0,     // eday of epoc
		39.48168677, // semi major axis (au)
		0.24880766,  // eccentricity
		17.14175,    // inclination (deg)
		110.30347,   // longitude of the ascending node (deg)
		224.06676,   // longitude of perihelion (deg)
		238.92881,   // mean longitude (deg)
		-0.00076912, // (au/julian century)
		0.00006465,  // (e/julian century)
		11.07,       // (arcsec/julian century)
		-37.33,      // (arcsec/julian century)
		-132.25,     // (arcsec/julian century)
		522747.90,   // (arcsec/julian century)
	}
	solstr = []string{"Fall equinox", "Winter solstice", "Spring equinox", "Summer solstice"}
)

func usage() {
	fmt.Fprintf(os.Stderr, "usage: astro [-jpokm] [-c nperiod] [-C tperiod] [-d date] [-e obj1 obj2] [-l nlat wlong elev] [-t ΔT]\n")
	os.Exit(2)
}

func main() {
	log.SetPrefix("astro: ")
	log.SetFlags(0)
	flag.Usage = usage
	flag.Parse()
	if flag.NArg() != 0 {
		usage()
	}
	t := time.Now().UTC()
	var err error
	if *startDate != "" {
		t, err = time.Parse(time.RFC3339, *startDate)
		if err != nil {
			log.Fatal(err)
		}
	}
	day = timeToJulian(&t)
	if *julian {
		fmt.Printf("Julian date: %.4f\n", day)
	}
	ΔT = deltaT(day)
	ostar = obj2{name: "Star", f: star}
	var eobjs [2]obj2
	if *eclipse != "" {
		if _, err := fmt.Sscanf(*eclipse, "%s %s", &eobjs[0].name, &eobjs[1].name); err != nil {
			log.Fatal("failed to parse eclipse objects")
		}
		for i, e := range eobjs {
			j := slices.IndexFunc(objs, func(o *obj2) bool {
				return e.name == o.name || e.name == o.fname
			})
			if j < 0 {
				log.Fatal("failed to parse eclipse objects")
			}
			eobjs[i] = *objs[j]
		}
	}
	root = os.Getenv("PLAN9")
	if root == "" {
		root = "/usr/local/plan9"
	}
	// Murray Hill, NJ
	nlat = (40 + 41.06/60) * radian
	wlong = (74 + 23.98/60) * radian
	elev = 150 * metersToFeet
	if *initLoc != "" {
		if err = parseLocation(*initLoc); err != nil {
			log.Fatal("failed to parse location")
		}
	} else {
		data, err := os.ReadFile(filepath.Join(root, "sky", "here"))
		if err == nil {
			_ = parseLocation(string(data))
		}
	}
	glat = nlat - (692.74*radsec)*math.Sin(2*nlat) + (1.16*radsec)*math.Sin(4*nlat)
	erad = .99832707e0 + .00167644e0*math.Cos(2*nlat) - .352e-5*math.Cos(4*nlat) + .001e-5*math.Cos(6*nlat) + .1568e-6*elev
	for range *nperiods {
		d := day
		fmt.Print(julianToTime(d))
		if *pflag || *eclipse != "" {
			pstime(d)
		}
		fmt.Println()
		for i := range objs[0].point {
			setime(d)
			for j := range objs {
				objs[j].f()
				setobj(&objs[j].point[i])
				if *pflag {
					if *mflag && objs[j].name == "Comet" {
						continue
					}
					output(objs[j].name, objs[j].point[i])
				}
			}
			if *eclipse != "" {
				d = dist(eobjs[0].point[i], eobjs[1].point[i])
				fmt.Printf("dist %s to %s = %.4f\n", eobjs[0].name, eobjs[1].name, d)
			}
			if *pflag || *eclipse != "" {
				break
			}
			d += deld
		}
		if !*pflag || *eclipse == "" {
			if err := search(); err != nil {
				log.Fatal(err)
			}
		}
		day += *p
	}
}

var t1899 = time.Date(1899, 12, 31, 12, 0, 0, 0, time.UTC)

func timeToJulian(t *time.Time) float64 {
	if *local {
		*t = t.Local()
	}
	return float64(t.Unix()-t1899.Unix()) / secondsPerDay
}

func julianToTime(jd float64) time.Time {
	s := jd * secondsPerDay
	t := t1899.Add(time.Duration(s) * time.Second)
	if *local {
		return t.Local()
	}
	return t
}

func deltaT(jDay float64) float64 {
	if *dt != 0 {
		return *dt
	}
	return jDay * meanSolarDaySeconds
}

func parseLocation(s string) error {
	if _, err := fmt.Sscanf(s, "%f %f %f", &nlat, &awlong, &elev); err != nil {
		return err
	}
	nlat *= radian
	awlong *= radian
	elev *= metersToFeet
	return nil
}

func fsun() {
	beta = 0
	rad = 0
	lambda = 0
	motion = 0
	helio()
	geo()
	seday = eday
	salph = alpha
	sdelt = delta
	mag = lmb2
}

func sun() {
	var mven, merth, mmars, mjup, msat, dmoon, mmoon, gmoon,
		pturbb, pturbl, pturbr, lograd float64
	ecc = .01675104 - 4.180e-5*capt - 1.26e-7*capt2
	incl = 0
	node = 0
	argp = 281.220833 + .0000470684*eday + .000453*capt2 + .000003*capt3
	mrad = 1
	anom = 358.475845 + .9856002670*eday - .000150*capt2 - .000003*capt3
	motion = .9856473354
	dmoon = 350.737681 + 12.1907491914*eday - .001436*capt2
	gmoon = 11.250889 + 13.2293504490*eday - .003212*capt2
	mmoon = 296.104608 + 13.0649924465*eday + 9.192e-3*capt2
	mven = 212.448 + 1.602121635*eday
	merth = 358.476 + 0.985600267*eday
	mmars = 319.590 + .524024095*eday
	mjup = 225.269 + .083082362*eday
	msat = 175.593 + .033450794*eday
	dmoon = math.Mod(dmoon, 360.) * radian
	gmoon = math.Mod(gmoon, 360.) * radian
	mmoon = math.Mod(mmoon, 360.) * radian
	mven *= radian
	merth *= radian
	mmars *= radian
	mjup *= radian
	msat *= radian
	anom += cosadd(sunf[0], sunc[0], []float64{mmars, merth, mven, mjup}) / 3600.
	anom += sinadd(sunf[1], sunc[1], []float64{mmars, merth, mven, mjup, .07884 * capt}) / 3600.
	incl *= radian
	node *= radian
	argp *= radian
	anom = math.Mod(anom, 360.) * radian
	// Computation of elliptic orbit.
	lambda = anom + argp
	pturbl = (6910.057-17.240*capt-0.052*capt2)*math.Sin(anom) +
		(72.338-0.361*capt)*math.Sin(2.*anom) +
		(1.054-0.001*capt)*math.Sin(3.*anom) + 0.018*math.Sin(4.*anom)
	lambda += pturbl * radsec
	beta = 0.
	lograd = (30.57e-6 - 0.15e-6*capt) -
		(7274.12e-6-18.14e-6*capt-0.05e-6*capt2)*math.Cos(anom) -
		(91.38e-6-0.46e-6*capt)*math.Cos(2.*anom) -
		(1.45e-6-0.01e-6*capt)*math.Cos(3.*anom) -
		0.02e-6*math.Cos(4.*anom)
	pturbl = cosadd(sunf[2], sunc[2], []float64{mmars, merth, mven, mjup, msat})
	pturbl += sinadd(sunf[3], sunc[3], []float64{dmoon, mmoon, merth}) + .9
	pturbl *= radsec
	pturbb = cosadd(sunf[4], sunc[4], []float64{merth, mven, mjup})
	pturbb += sinadd(sunf[5], sunc[5], []float64{gmoon, mmoon, dmoon})
	pturbb *= radsec
	pturbr = cosadd(sunf[6], sunc[6], []float64{mmars, merth, mven, mjup, msat})
	pturbr += cosadd(sunf[7], sunc[7], []float64{dmoon, mmoon, merth})
	lambda += pturbl
	if lambda > pipi {
		lambda -= pipi
	}
	beta += pturbb
	lograd = (lograd + pturbr) * 2.30258509
	rad = 1 + lograd*(1+lograd*(.5+lograd/6))
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 961.182
	if *oflag {
		semi = 959.63
	}
	mag = -26.5
}

// helio converts from ecliptic heliocentric coordinates referred to the mean
// equinox of date to equatorial geocentric coordinates referred to the true
// equator and equinox.
func helio() {
	var xmp, ymp, zmp, beta2 float64
	// Compute geocentric distance of object and light-time correction.
	xmp = rad * math.Cos(beta) * math.Cos(lambda)
	ymp = rad * math.Cos(beta) * math.Sin(lambda)
	zmp = rad * math.Sin(beta)
	rp = math.Sqrt((xmp+xms)*(xmp+xms) + (ymp+yms)*(ymp+yms) + (zmp+zms)*(zmp+zms))
	lmb2 = lambda - .0057756e0*rp*motion
	xmp = rad * math.Cos(beta) * math.Cos(lmb2)
	ymp = rad * math.Cos(beta) * math.Sin(lmb2)
	zmp = rad * math.Sin(beta)
	// Compute annual parallax from the position of the sun.
	xmp += xms
	ymp += yms
	zmp += zms
	rp = math.Sqrt(xmp*xmp + ymp*ymp + zmp*zmp)
	// Compute annual aberration from the orbital velocity of the earth (by an incorrect method).
	xmp -= xdot * rp
	ymp -= ydot * rp
	zmp -= zdot * rp
	// Perform the nutation and so convert from the mean equator and equinox to the true.
	lmb2 = math.Atan2(ymp, xmp)
	beta2 = math.Atan2(zmp, math.Sqrt(xmp*xmp+ymp*ymp))
	lmb2 += phi
	// Change to equatorial coordinates.
	xmp = rp * math.Cos(lmb2) * math.Cos(beta2)
	ymp = rp * (math.Sin(lmb2)*math.Cos(beta2)*math.Cos(tobliq) - math.Sin(tobliq)*math.Sin(beta2))
	zmp = rp * (math.Sin(lmb2)*math.Cos(beta2)*math.Sin(tobliq) + math.Cos(tobliq)*math.Sin(beta2))
	alpha = math.Atan2(ymp, xmp)
	delta = math.Atan2(zmp, math.Sqrt(xmp*xmp+ymp*ymp))
	hp = 8.794e0 * radsec / rp
	semi /= rp
	if rad > 0 && rad < 2.e5 {
		mag += 2.17 * math.Log(rad*rp)
	}
}

// geo converts geocentric equatorial coordinates to topocentric equatorial and
// topocentric horizon coordinates. All are (usually) referred to the true equator.
func geo() {
	var sel, saz, caz, f, sa, ca, sd float64
	// Convert to local hour angle and declination
	lha = gst - alpha - wlong
	decl = delta
	// Compute diurnal parallax (requires geocentric latitude)
	sa = math.Cos(decl) * math.Sin(lha)
	ca = math.Cos(decl)*math.Cos(lha) - erad*math.Cos(glat)*math.Sin(hp)
	sd = math.Sin(decl) - erad*math.Sin(glat)*math.Sin(hp)
	lha = math.Atan2(sa, ca)
	decl2 = math.Atan2(sd, math.Sqrt(sa*sa+ca*ca))
	f = math.Sqrt(sa*sa + ca*ca + sd*sd)
	semi2 = semi / f
	ra = gst - lha - wlong
	ra = pinorm(ra)
	// Convert to horizon coordinates
	sel = math.Sin(nlat)*math.Sin(decl2) + math.Cos(nlat)*math.Cos(decl2)*math.Cos(lha)
	el = math.Atan2(sel, pyth(sel))
	saz = math.Sin(lha) * math.Cos(decl2)
	caz = math.Cos(nlat)*math.Sin(decl2) - math.Sin(nlat)*math.Cos(decl2)*math.Cos(lha)
	az = math.Pi + math.Atan2(saz, -caz)
	az /= radian
	el /= radian
}

func pinorm(a float64) float64 {
	for a < -math.Pi {
		a += pipi
	}
	for a > math.Pi {
		a -= pipi
	}
	return a
}

func pyth(x float64) float64 {
	x *= x
	if x > 1 {
		x = 1
	}
	return math.Sqrt(1 - x)
}

var k1, k2, k3, k4, mnom, msun, noded, dmoon float64

func moon() {
	var dlong, lsun, psun, eccm, eccs, chp, cpe, v0, t0, m0, j0,
		arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10,
		dgamma, k5, k6, lterms, sterms, cterms, nterms, pterms, spterms,
		gamma1, gamma2, gamma3, arglat, xmp, ymp, zmp, obl2 float64
	// The fundamental elements - all referred to the epoch of Jan 0.5, 1900 and
	// to the mean equinox of date.
	dlong = 270.434164 + 13.1763965268*eday - .001133*capt2 + 2.e-6*capt3
	argp = 334.329556 + .1114040803*eday - .010325*capt2 - 12.e-6*capt3
	node = 259.183275 - .0529539222*eday + .002078*capt2 + 2.e-6*capt3
	lsun = 279.696678 + .9856473354*eday + .000303*capt2
	psun = 281.220833 + .0000470684*eday + .000453*capt2 + 3.e-6*capt3
	dlong = math.Mod(dlong, 360.)
	argp = math.Mod(argp, 360.)
	node = math.Remainder(node, 360.)
	lsun = math.Mod(lsun, 360.)
	psun = math.Mod(psun, 360.)
	eccm = 22639.550
	eccs = .01675104 - .00004180*capt
	incl = 18461.400
	cpe = 124.986
	chp = 3422.451
	// Some subsidiary elements - they are all longitudes and they are referred
	// to the epoch 1/0.5 1900 and to the fixed mean equinox of 1850.0.
	v0 = 342.069128 + 1.6021304820*eday
	t0 = 98.998753 + 0.9856091138*eday
	m0 = 293.049675 + 0.5240329445*eday
	j0 = 237.352319 + 0.0830912295*eday
	// The following are periodic corrections to the fundamental elements and constants.
	arg1 = 41.1 + 20.2*(capt+.5)
	arg2 = dlong - argp + 33. + 3.*t0 - 10.*v0 - 2.6*(capt+.5)
	arg3 = dlong - argp + 151.1 + 16.*t0 - 18.*v0 - (capt + .5) // The "Great Venus Inequality".
	arg4 = node
	arg5 = node + 276.2 - 2.3*(capt+.5)
	arg6 = 313.9 + 13.*t0 - 8.*v0
	arg7 = dlong - argp + 112.0 + 29.*t0 - 26.*v0
	arg8 = dlong + argp - 2.*lsun + 273. + 21.*t0 - 20.*v0
	arg9 = node + 290.1 - 0.9*(capt+.5)
	arg10 = 115. + 38.5*(capt+.5)
	arg1 *= radian
	arg2 *= radian
	arg3 *= radian
	arg4 *= radian
	arg5 *= radian
	arg6 *= radian
	arg7 *= radian
	arg8 *= radian
	arg9 *= radian
	arg10 *= radian
	dlong += (0.84*math.Sin(arg1) + 0.31*math.Sin(arg2) + 14.27*math.Sin(arg3) + 7.261*math.Sin(arg4) +
		0.282*math.Sin(arg5) + 0.237*math.Sin(arg6) + 0.108*math.Sin(arg7) + 0.126*math.Sin(arg8)) / 3600.
	argp += (-2.10*math.Sin(arg1) - 0.118*math.Sin(arg3) - 2.076*math.Sin(arg4) - 0.840*math.Sin(arg5) - 0.593*math.Sin(arg6)) / 3600.
	node += (0.63*math.Sin(arg1) + 0.17*math.Sin(arg3) + 95.96*math.Sin(arg4) + 15.58*math.Sin(arg5) + 1.86*math.Sin(arg9)) / 3600.
	t0 += (-6.40*math.Sin(arg1) - 1.89*math.Sin(arg6)) / 3600.
	psun += (6.40*math.Sin(arg1) + 1.89*math.Sin(arg6)) / 3600.
	dgamma = -4.318*math.Cos(arg4) - 0.698*math.Cos(arg5) - 0.083*math.Cos(arg9)
	j0 += 0.33 * math.Sin(arg10)
	// The following factors account for the fact that the eccentricity, solar eccentricity,
	// inclination and parallax used by Brown to make up his coefficients are
	// both wrong and out of date.  Brown did the same thing in a different way.
	k1 = eccm / 22639.500
	k2 = eccs / .01675104
	k3 = 1. + 2.708e-6 + .000108008*dgamma
	k4 = cpe / 125.154
	k5 = chp / 3422.700
	// The principal arguments that are used to compute perturbations are the following
	// differences of the fundamental elements.
	mnom = dlong - argp
	msun = lsun - psun
	noded = dlong - node
	dmoon = dlong - lsun
	// Solar terms in longitude.
	var i int
	for ; moontab[i].f != 0.0; i++ {
		lterms += sinx(moontab[i].f, moontab[i].c[0], moontab[i].c[1], moontab[i].c[2], moontab[i].c[3], 0.0)
	}
	i++
	fmt.Printf("first lterms=%f, i=%d\n", lterms, i)
	// Planetary terms in longitude.
	lterms += sinx(0.822, 0, 0, 0, 0, t0-v0)
	lterms += sinx(0.307, 0, 0, 0, 0, 2.*t0-2.*v0+179.8)
	lterms += sinx(0.348, 0, 0, 0, 0, 3.*t0-2.*v0+272.9)
	lterms += sinx(0.176, 0, 0, 0, 0, 4.*t0-3.*v0+271.7)
	lterms += sinx(0.092, 0, 0, 0, 0, 5.*t0-3.*v0+199.)
	lterms += sinx(0.129, 1, 0, 0, 0, -t0+v0+180.)
	lterms += sinx(0.152, 1, 0, 0, 0, t0-v0)
	lterms += sinx(0.127, 1, 0, 0, 0, 3.*t0-3.*v0+180.)
	lterms += sinx(0.099, 0, 0, 0, 2, t0-v0)
	lterms += sinx(0.136, 0, 0, 0, 2, 2.*t0-2.*v0+179.5)
	lterms += sinx(0.083, -1, 0, 0, 2, -4.*t0+4.*v0+180.)
	lterms += sinx(0.662, -1, 0, 0, 2, -3.*t0+3.*v0+180.0)
	lterms += sinx(0.137, -1, 0, 0, 2, -2.*t0+2.*v0)
	lterms += sinx(0.133, -1, 0, 0, 2, t0-v0)
	lterms += sinx(0.157, -1, 0, 0, 2, 2.*t0-2.*v0+179.6)
	lterms += sinx(0.079, -1, 0, 0, 2, -8.*t0+6.*v0+162.6)
	lterms += sinx(0.073, 2, 0, 0, -2, 3.*t0-3.*v0+180.)
	lterms += sinx(0.643, 0, 0, 0, 0, -t0+j0+178.8)
	lterms += sinx(0.187, 0, 0, 0, 0, -2.*t0+2.*j0+359.6)
	lterms += sinx(0.087, 0, 0, 0, 0, j0+289.9)
	lterms += sinx(0.165, 0, 0, 0, 0, -t0+2.*j0+241.5)
	lterms += sinx(0.144, 1, 0, 0, 0, t0-j0+1.0)
	lterms += sinx(0.158, 1, 0, 0, 0, -t0+j0+179.0)
	lterms += sinx(0.190, 1, 0, 0, 0, -2.*t0+2.*j0+180.0)
	lterms += sinx(0.096, 1, 0, 0, 0, -2.*t0+3.*j0+352.5)
	lterms += sinx(0.070, 0, 0, 0, 2, 2.*t0-2.*j0+180.)
	lterms += sinx(0.167, 0, 0, 0, 2, -t0+j0+178.5)
	lterms += sinx(0.085, 0, 0, 0, 2, -2.*t0+2.*j0+359.2)
	lterms += sinx(1.137, -1, 0, 0, 2, 2.*t0-2.*j0+180.3)
	lterms += sinx(0.211, -1, 0, 0, 2, -t0+j0+178.4)
	lterms += sinx(0.089, -1, 0, 0, 2, -2.*t0+2.*j0+359.2)
	lterms += sinx(0.436, -1, 0, 0, 2, 2.*t0-3.*j0+7.5)
	lterms += sinx(0.240, 2, 0, 0, -2, -2.*t0+2.*j0+179.9)
	lterms += sinx(0.284, 2, 0, 0, -2, -2.*t0+3.*j0+172.5)
	lterms += sinx(0.195, 0, 0, 0, 0, -2.*t0+2.*m0+180.2)
	lterms += sinx(0.327, 0, 0, 0, 0, -t0+2.*m0+224.4)
	lterms += sinx(0.093, 0, 0, 0, 0, -2.*t0+4.*m0+244.8)
	lterms += sinx(0.073, 1, 0, 0, 0, -t0+2.*m0+223.3)
	lterms += sinx(0.074, 1, 0, 0, 0, t0-2.*m0+306.3)
	lterms += sinx(0.189, 0, 0, 0, 0, node+180.)
	fmt.Printf("lterms=%f\n", lterms)
	// Solar terms in latitude.
	sterms = 0
	for {
		if moontab[i].f == 0 {
			break
		}
		sterms += sinx(moontab[i].f, moontab[i].c[0], moontab[i].c[1], moontab[i].c[2], moontab[i].c[3], 0)
		i++
	}
	i++
	cterms = 0
	for {
		if moontab[i].f == 0 {
			break
		}
		cterms += cosx(moontab[i].f, moontab[i].c[0], moontab[i].c[1], moontab[i].c[2], moontab[i].c[3], 0)
		i++
	}
	i++
	nterms = 0
	for {
		if moontab[i].f == 0 {
			break
		}
		nterms += sinx(moontab[i].f, moontab[i].c[0], moontab[i].c[1], moontab[i].c[2], moontab[i].c[3], 0)
		i++
	}
	i++
	// Planetary terms in latitude.
	pterms = sinx(0.215, 0, 0, 0, 0, dlong)
	// Solar terms in parallax.
	spterms = 3422.700
	for {
		if moontab[i].f == 0 {
			break
		}
		spterms += cosx(moontab[i].f, moontab[i].c[0], moontab[i].c[1], moontab[i].c[2], moontab[i].c[3], 0)
		i++
	}
	// Planetary terms in parallax.
	// spterms = spterms
	// Computation of longitude.
	lambda = (dlong + lterms/3600.) * radian
	// Computation of latitude.
	arglat = (noded + sterms/3600.) * radian
	gamma1 = 18519.700 * k3
	gamma2 = -6.241 * k3 * k3 * k3
	gamma3 = 0.004 * k3 * k3 * k3 * k3 * k3
	k6 = (gamma1 + cterms) / gamma1
	beta = k6*(gamma1*math.Sin(arglat)+gamma2*math.Sin(3.*arglat)+gamma3*math.Sin(5.*arglat)+nterms) + pterms
	if *oflag {
		beta -= 0.6
	}
	beta *= radsec
	// Computation of parallax.
	spterms = k5 * spterms * radsec
	hp = spterms + (spterms*spterms*spterms)/6.
	rad = hp / radsec
	rp = 1.
	semi = .0799 + .272453*(hp/radsec)
	if dmoon < 0. {
		dmoon += 360.
	}
	mag = dmoon / 360.
	// Change to equatorial coordinates.
	lambda += phi
	obl2 = obliq + eps
	xmp = rp * math.Cos(lambda) * math.Cos(beta)
	ymp = rp * (math.Sin(lambda)*math.Cos(beta)*math.Cos(obl2) - math.Sin(obl2)*math.Sin(beta))
	zmp = rp * (math.Sin(lambda)*math.Cos(beta)*math.Sin(obl2) + math.Cos(obl2)*math.Sin(beta))
	alpha = math.Atan2(ymp, xmp)
	delta = math.Atan2(zmp, math.Sqrt(xmp*xmp+ymp*ymp))
	meday = eday
	mhp = hp
	geo()
}

func sinx(coef float64, i, j, k, m int, angle float64) float64 {
	return trigx(coef, i, j, k, m, angle, math.Sin)
}

func cosx(coef float64, i, j, k, m int, angle float64) float64 {
	return trigx(coef, i, j, k, m, angle, math.Cos)
}

func trigx(coef float64, i, j, k, m int, angle float64, f func(float64) float64) float64 {
	x := coef * f((float64(i)*mnom+float64(j)*msun+float64(k)*noded+float64(m)*dmoon+angle)*radian)
	for i := abs(i); i > 0; i-- {
		x *= k1
	}
	for j := abs(j); j > 0; j-- {
		x *= k2
	}
	for k := abs(k); k > 0; k-- {
		x *= k3
	}
	if m&1 != 0 {
		x *= k4
	}
	return x
}

func abs(x int) int {
	if x < 0 {
		x = -x
	}
	return x
}

func shad() {
	if seday != eday {
		fsun()
	}
	if meday != eday {
		moon()
	}
	alpha = math.Mod(salph+math.Pi, pipi)
	delta = -sdelt
	hp = mhp
	semi = 1.0183*mhp/radsec - 969.85/srad
	geo()
}

func merc() {
	var pturbl, pturbr, lograd, dele, enom, vnom, nd, sl, q0, v0, t0, j0, s0,
		lsun, elong, ci, dlong float64
	ecc = .20561421 + .00002046*capt - 0.03e-6*capt2
	incl = 7.0028806 + .0018608*capt - 18.3e-6*capt2
	node = 47.145944 + 1.185208*capt + .0001739*capt2
	argp = 75.899697 + 1.555490*capt + .0002947*capt2
	mrad = .3870986
	anom = 102.279381 + 4.0923344364*eday + 6.7e-6*capt2
	motion = 4.0923770233
	q0 = 102.28 + 4.092334429*eday
	v0 = 212.536 + 1.602126105*eday
	t0 = -1.45 + .985604737*eday
	j0 = 225.36 + .083086735*eday
	s0 = 175.68 + .033455441*eday
	q0 *= radian
	v0 *= radian
	t0 *= radian
	j0 *= radian
	s0 *= radian
	incl *= radian
	node *= radian
	argp *= radian
	anom = math.Mod(anom, 360.) * radian
	enom = anom + ecc*math.Sin(anom)
	for {
		dele = (anom - enom + ecc*math.Sin(enom)) / (1. - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom = 2. * math.Atan2(math.Sqrt((1.+ecc)/(1.-ecc))*math.Sin(enom/2.), math.Cos(enom/2.))
	rad = mrad * (1. - ecc*math.Cos(enom))
	pturbl = cosadd(mercf[0], mercc[0], []float64{q0, -v0})
	pturbl += cosadd(mercf[1], mercc[1], []float64{q0, -t0})
	pturbl += cosadd(mercf[2], mercc[2], []float64{q0, -j0})
	pturbl += cosadd(mercf[3], mercc[3], []float64{q0, -s0})
	pturbr = cosadd(mercf[4], mercc[4], []float64{q0, -v0})
	pturbr += cosadd(mercf[5], mercc[5], []float64{q0, -t0})
	pturbr += cosadd(mercf[6], mercc[6], []float64{q0, -j0})
	// Reduce to the ecliptic.
	lambda = vnom + argp + pturbl*radsec
	nd = lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl = math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, pyth(sl))
	lograd = pturbr * 2.30258509
	rad *= 1. + lograd
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 3.34
	lsun = 99.696678 + 0.9856473354*eday
	lsun *= radian
	elong = lambda - lsun
	ci = (rad - math.Cos(elong)) / math.Sqrt(1.+rad*rad-2.*rad*math.Cos(elong))
	dlong = math.Atan2(pyth(ci), ci) / radian
	mag = -.003 + .01815*dlong + .0001023*dlong*dlong
	helio()
	geo()
}

func cosadd(caf []float64, cac []int, coefs []float64) float64 {
	return trigadd(math.Cos, caf, cac, coefs)
}

func sinadd(caf []float64, cac []int, coefs []float64) float64 {
	return trigadd(math.Sin, caf, cac, coefs)
}

func trigadd(f func(float64) float64, caf []float64, cac []int, coefs []float64) float64 {
	var sum float64
	var j int
	for i := 0; i < len(caf); i += 2 {
		a1, a2 := caf[i], caf[i+1]
		for _, c := range coefs {
			a2 += float64(cac[j]) * c
			j++
		}
		sum += a1 * f(a2)
	}
	return sum
}

func venus() {
	var pturbl, pturbb, pturbr, lograd, dele, enom, vnom, nd, sl, v0, t0, m0, j0, s0,
		lsun, elong, ci, dlong float64
	// Mean orbital elements.
	ecc = .00682069 - .00004774*capt + 0.091e-6*capt2
	incl = 3.393631 + .0010058*capt - 0.97e-6*capt2
	node = 75.779647 + .89985*capt + .00041*capt2
	argp = 130.163833 + 1.408036*capt - .0009763*capt2
	mrad = .7233316
	anom = 212.603219 + 1.6021301540*eday + .00128605*capt2
	motion = 1.6021687039
	// Mean anomalies of perturbing planets.
	v0 = 212.60 + 1.602130154*eday
	t0 = 358.63 + .985608747*eday
	m0 = 319.74 + 0.524032490*eday
	j0 = 225.43 + .083090842*eday
	s0 = 175.8 + .033459258*eday
	v0 *= radian
	t0 *= radian
	m0 *= radian
	j0 *= radian
	s0 *= radian
	incl *= radian
	node *= radian
	argp *= radian
	anom = math.Mod(anom, 360.) * radian
	// Computation of long period terms affecting the mean anomaly.
	anom += (2.761-0.022*capt)*radsec*math.Sin(13.*t0-8.*v0+43.83*radian+4.52*radian*capt) +
		0.268*radsec*math.Cos(4.*m0-7.*t0+3.*v0) +
		0.019*radsec*math.Sin(4.*m0-7.*t0+3.*v0) -
		0.208*radsec*math.Sin(s0+1.4*radian*capt)
	// Computation of elliptic orbit.
	enom = anom + ecc*math.Sin(anom)
	for {
		dele = (anom - enom + ecc*math.Sin(enom)) / (1. - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom = 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	lambda = vnom + argp
	// Perturbations in longitude.
	pturbl = cosadd(venf[0], venc[0], []float64{v0, t0, m0, j0})
	pturbl *= radsec
	// Perturbations in latitude.
	pturbb = cosadd(venf[1], venc[1], []float64{v0, t0, j0})
	pturbb *= radsec
	// Perturbations in log radius vector.
	pturbr = cosadd(venf[2], venc[2], []float64{v0, t0, m0, j0})
	// Reduction to the ecliptic.
	lambda += pturbl
	nd = lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl = math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, pyth(sl)) + pturbb
	lograd = pturbr * 2.30258509
	rad *= 1 + lograd
	motion *= radian * mrad * mrad / (rad * rad)
	// Computation of magnitude.
	lsun = 99.696678 + 0.9856473354*eday
	lsun *= radian
	elong = lambda - lsun
	ci = (rad - math.Cos(elong)) / math.Sqrt(1+rad*rad-2*rad*math.Cos(elong))
	dlong = math.Atan2(pyth(ci), ci) / radian
	mag = -4 + .01322*dlong + .0000004247*dlong*dlong*dlong
	semi = 8.41
	helio()
	geo()
}

func mars() {
	var pturbl, pturbb, pturbr, lograd, dele, enom, vnom, nd, sl,
		lsun, elong, ci, dlong float64
	ecc = .09331290 + .000092064*capt
	incl = 1.850333 - 6.75e-4*capt
	node = 48.786442 + .770992*capt
	argp = 334.218203 + 1.840758*capt + 1.30e-4*capt2
	mrad = 1.5236915
	anom = 319.529425 + .5240207666*eday + 1.808e-4*capt2
	motion = 0.5240711638
	incl = incl * radian
	node = node * radian
	argp = argp * radian
	anom = math.Mod(anom, 360.) * radian
	enom = anom + ecc*math.Sin(anom)
	for {
		dele = (anom - enom + ecc*math.Sin(enom)) / (1. - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom = 2. * math.Atan2(math.Sqrt((1.+ecc)/(1.-ecc))*math.Sin(enom/2.), math.Cos(enom/2.))
	rad = mrad * (1. - ecc*math.Cos(enom))
	lambda = vnom + argp
	pturbl = 0.
	lambda = lambda + pturbl*radsec
	pturbb = 0.
	pturbr = 0.
	// Reduce to the ecliptic.
	nd = lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl = math.Sin(incl)*math.Sin(nd) + pturbb*radsec
	beta = math.Atan2(sl, pyth(sl))
	lograd = pturbr * 2.30258509
	rad *= 1. + lograd
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 4.68
	lsun = 99.696678 + 0.9856473354*eday
	lsun *= radian
	elong = lambda - lsun
	ci = (rad - math.Cos(elong)) / math.Sqrt(1.+rad*rad-2.*rad*math.Cos(elong))
	dlong = math.Atan2(pyth(ci), ci) / radian
	mag = -1.30 + .01486*dlong
	helio()
	geo()
}

func jup() {
	var pturbl, pturbb, pturbr, lograd, dele, enom, vnom, nd, sl float64
	ecc = .0483376 + 163.e-6*capt
	incl = 1.308660 - .0055*capt
	node = 99.43785 + 1.011*capt
	argp = 12.71165 + 1.611*capt
	mrad = 5.202803
	anom = 225.22165 + .0830912*eday - .0484*capt
	motion = 299.1284 / 3600.
	incl *= radian
	node *= radian
	argp *= radian
	anom = math.Mod(anom, 360.) * radian
	enom = anom + ecc*math.Sin(anom)
	for {
		dele = (anom - enom + ecc*math.Sin(enom)) / (1. - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom = 2. * math.Atan2(math.Sqrt((1.+ecc)/(1.-ecc))*math.Sin(enom/2.), math.Cos(enom/2.))
	rad = mrad * (1. - ecc*math.Cos(enom))
	lambda = vnom + argp
	pturbl = 0.
	lambda += pturbl * radsec
	pturbb = 0.
	pturbr = 0.
	// Reduce to the ecliptic.
	nd = lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl = math.Sin(incl)*math.Sin(nd) + pturbb*radsec
	beta = math.Atan2(sl, pyth(sl))
	lograd = pturbr * 2.30258509
	rad *= 1. + lograd
	lambda += 555. * radsec
	beta -= 51. * radsec
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 98.47
	mag = -8.93
	helio()
	geo()
}

func sat() {
	var pturbl, pturbb, pturbr, lograd, dele, enom, vnom, nd, sl,
		capj, capn, eye, comg, omg, sb, su, cu, u, b, up, sd, ca, sa float64
	ecc = .0558900 - .000347*capt
	incl = 2.49256 - .0044*capt
	node = 112.78364 + .87306*capt
	argp = 91.08897 + 1.95917*capt
	mrad = 9.538843
	anom = 175.47630 + .03345972*eday - .56527*capt
	motion = 120.4550 / 3600.
	incl *= radian
	node *= radian
	argp *= radian
	anom = math.Mod(anom, 360.) * radian
	enom = anom + ecc*math.Sin(anom)
	for {
		dele = (anom - enom + ecc*math.Sin(enom)) / (1. - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom = 2. * math.Atan2(math.Sqrt((1.+ecc)/(1.-ecc))*math.Sin(enom/2.), math.Cos(enom/2.))
	rad = mrad * (1. - ecc*math.Cos(enom))
	lambda = vnom + argp
	pturbl = 0.
	lambda += pturbl * radsec
	pturbb = 0.
	pturbr = 0.
	// Reduce to the ecliptic.
	nd = lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl = math.Sin(incl)*math.Sin(nd) + pturbb*radsec
	beta = math.Atan2(sl, pyth(sl))
	lograd = pturbr * 2.30258509
	rad *= 1. + lograd
	lambda -= 1185. * radsec
	beta -= 51. * radsec
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 83.33
	// Computation of magnitude; first, find the geocentric equatorial coordinates of Saturn.
	sd = rad * (math.Cos(beta)*math.Sin(lambda)*math.Sin(obliq) + math.Sin(beta)*math.Cos(obliq))
	sa = rad * (math.Cos(beta)*math.Sin(lambda)*math.Cos(obliq) - math.Sin(beta)*math.Sin(obliq))
	ca = rad * math.Cos(beta) * math.Cos(lambda)
	sd += zms
	sa += yms
	ca += xms
	alpha = math.Atan2(sa, ca)
	delta = math.Atan2(sd, math.Sqrt(sa*sa+ca*ca))
	// Here are the necessary elements of Saturn's rings cf. Exp. Supp. p. 363ff.
	capj = 6.9056 - 0.4322*capt
	capn = 126.3615 + 3.9894*capt + 0.2403*capt2
	eye = 28.0743 - 0.0128*capt
	comg = 168.1179 + 1.3936*capt
	omg = 42.9236 - 2.7390*capt - 0.2344*capt2
	capj *= radian
	capn *= radian
	eye *= radian
	comg *= radian
	omg *= radian
	// Now find saturnicentric ring-plane coords of the earth.
	sb = math.Sin(capj)*math.Cos(delta)*math.Sin(alpha-capn) - math.Cos(capj)*math.Sin(delta)
	su = math.Cos(capj)*math.Cos(delta)*math.Sin(alpha-capn) + math.Sin(capj)*math.Sin(delta)
	cu = math.Cos(delta) * math.Cos(alpha-capn)
	u = math.Atan2(su, cu)
	b = math.Atan2(sb, math.Sqrt(su*su+cu*cu))
	// And then the saturnicentric ring-plane coords of the sun.
	su = math.Sin(eye)*math.Sin(beta) + math.Cos(eye)*math.Cos(beta)*math.Sin(lambda-comg)
	cu = math.Cos(beta) * math.Cos(lambda-comg)
	up = math.Atan2(su, cu)
	// At last, the magnitude.
	sb = math.Sin(b)
	mag = -8.68 + 2.52*math.Abs(up+omg-u) - 2.60*math.Abs(sb) + 1.25*(sb*sb)
	helio()
	geo()
}

func uran() {
	var pturbl, pturbb, pturbr, lograd, dele, enom, vnom, nd, sl,
		capj, capn, eye, comg, omg, sb, su, cu, u, b, up, sd, ca, sa, cy float64
	cy = (eday - elemUran[0]) / 36525. // Per julian century
	mrad = elemUran[1] + elemUran[1+6]*cy
	ecc = elemUran[2] + elemUran[2+6]*cy
	cy = cy / 3600 // arcsec/deg per julian century
	incl = elemUran[3] + elemUran[3+6]*cy
	node = elemUran[4] + elemUran[4+6]*cy
	argp = elemUran[5] + elemUran[5+6]*cy
	anom = elemUran[6] + elemUran[6+6]*cy - argp
	motion = elemUran[6+6] / 36525. / 3600
	incl *= radian
	node *= radian
	argp *= radian
	anom = math.Mod(anom, 360.) * radian
	enom = anom + ecc*math.Sin(anom)
	for {
		dele = (anom - enom + ecc*math.Sin(enom)) / (1. - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom = 2. * math.Atan2(math.Sqrt((1.+ecc)/(1.-ecc))*math.Sin(enom/2.), math.Cos(enom/2.))
	rad = mrad * (1. - ecc*math.Cos(enom))
	lambda = vnom + argp
	pturbl = 0.
	lambda += pturbl * radsec
	pturbb = 0.
	pturbr = 0.
	// Reduce to the ecliptic.
	nd = lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl = math.Sin(incl)*math.Sin(nd) + pturbb*radsec
	beta = math.Atan2(sl, pyth(sl))
	lograd = pturbr * 2.30258509
	rad *= 1. + lograd
	lambda -= 1185. * radsec
	beta -= 51. * radsec
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 83.33
	// Computation of magnitude; first, find the geocentric equatorial coordinates of Saturn.
	sd = rad * (math.Cos(beta)*math.Sin(lambda)*math.Sin(obliq) + math.Sin(beta)*math.Cos(obliq))
	sa = rad * (math.Cos(beta)*math.Sin(lambda)*math.Cos(obliq) - math.Sin(beta)*math.Sin(obliq))
	ca = rad * math.Cos(beta) * math.Cos(lambda)
	sd += zms
	sa += yms
	ca += xms
	alpha = math.Atan2(sa, ca)
	delta = math.Atan2(sd, math.Sqrt(sa*sa+ca*ca))
	// Here are the necessary elements of Saturn's rings cf. Exp. Supp. p. 363ff.
	capj = 6.9056 - 0.4322*capt
	capn = 126.3615 + 3.9894*capt + 0.2403*capt2
	eye = 28.0743 - 0.0128*capt
	comg = 168.1179 + 1.3936*capt
	omg = 42.9236 - 2.7390*capt - 0.2344*capt2
	capj *= radian
	capn *= radian
	eye *= radian
	comg *= radian
	omg *= radian
	// Now find saturnicentric ring-plane coords of the earth.
	sb = math.Sin(capj)*math.Cos(delta)*math.Sin(alpha-capn) - math.Cos(capj)*math.Sin(delta)
	su = math.Cos(capj)*math.Cos(delta)*math.Sin(alpha-capn) + math.Sin(capj)*math.Sin(delta)
	cu = math.Cos(delta) * math.Cos(alpha-capn)
	u = math.Atan2(su, cu)
	b = math.Atan2(sb, math.Sqrt(su*su+cu*cu))
	// And then the saturnicentric ring-plane coords of the sun.
	su = math.Sin(eye)*math.Sin(beta) + math.Cos(eye)*math.Cos(beta)*math.Sin(lambda-comg)
	cu = math.Cos(beta) * math.Cos(lambda-comg)
	up = math.Atan2(su, cu)
	// At last, the magnitude.
	sb = math.Sin(b)
	mag = -8.68 + 2.52*math.Abs(up+omg-u) - 2.60*math.Abs(sb) + 1.25*(sb*sb)
	helio()
	geo()
}

func nept() {
	var pturbl, pturbb, pturbr, lograd, dele, enom, vnom, nd, sl,
		capj, capn, eye, comg, omg, sb, su, cu, u, b, up, sd, ca, sa, cy float64
	cy = (eday - elemNept[0]) / 36525. // Per julian century
	mrad = elemNept[1] + elemNept[1+6]*cy
	ecc = elemNept[2] + elemNept[2+6]*cy
	cy = cy / 3600 // arcsec/deg per julian century
	incl = elemNept[3] + elemNept[3+6]*cy
	node = elemNept[4] + elemNept[4+6]*cy
	argp = elemNept[5] + elemNept[5+6]*cy
	anom = elemNept[6] + elemNept[6+6]*cy - argp
	motion = elemNept[6+6] / 36525. / 3600
	incl *= radian
	node *= radian
	argp *= radian
	anom = math.Mod(anom, 360.) * radian
	enom = anom + ecc*math.Sin(anom)
	for {
		dele = (anom - enom + ecc*math.Sin(enom)) / (1. - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom = 2. * math.Atan2(math.Sqrt((1.+ecc)/(1.-ecc))*math.Sin(enom/2.), math.Cos(enom/2.))
	rad = mrad * (1. - ecc*math.Cos(enom))
	lambda = vnom + argp
	pturbl = 0.
	lambda += pturbl * radsec
	pturbb = 0.
	pturbr = 0.
	// Reduce to the ecliptic.
	nd = lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl = math.Sin(incl)*math.Sin(nd) + pturbb*radsec
	beta = math.Atan2(sl, pyth(sl))
	lograd = pturbr * 2.30258509
	rad *= 1. + lograd
	lambda -= 1185. * radsec
	beta -= 51. * radsec
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 83.33
	// Computation of magnitude; first, find the geocentric equatorial coordinates of Saturn.
	sd = rad * (math.Cos(beta)*math.Sin(lambda)*math.Sin(obliq) + math.Sin(beta)*math.Cos(obliq))
	sa = rad * (math.Cos(beta)*math.Sin(lambda)*math.Cos(obliq) - math.Sin(beta)*math.Sin(obliq))
	ca = rad * math.Cos(beta) * math.Cos(lambda)
	sd += zms
	sa += yms
	ca += xms
	alpha = math.Atan2(sa, ca)
	delta = math.Atan2(sd, math.Sqrt(sa*sa+ca*ca))
	// Here are the necessary elements of Saturn's rings cf. Exp. Supp. p. 363ff.
	capj = 6.9056 - 0.4322*capt
	capn = 126.3615 + 3.9894*capt + 0.2403*capt2
	eye = 28.0743 - 0.0128*capt
	comg = 168.1179 + 1.3936*capt
	omg = 42.9236 - 2.7390*capt - 0.2344*capt2
	capj *= radian
	capn *= radian
	eye *= radian
	comg *= radian
	omg *= radian
	// Now find saturnicentric ring-plane coords of the earth.
	sb = math.Sin(capj)*math.Cos(delta)*math.Sin(alpha-capn) - math.Cos(capj)*math.Sin(delta)
	su = math.Cos(capj)*math.Cos(delta)*math.Sin(alpha-capn) + math.Sin(capj)*math.Sin(delta)
	cu = math.Cos(delta) * math.Cos(alpha-capn)
	u = math.Atan2(su, cu)
	b = math.Atan2(sb, math.Sqrt(su*su+cu*cu))
	// And then the saturnicentric ring-plane coords of the sun.
	su = math.Sin(eye)*math.Sin(beta) + math.Cos(eye)*math.Cos(beta)*math.Sin(lambda-comg)
	cu = math.Cos(beta) * math.Cos(lambda-comg)
	up = math.Atan2(su, cu)
	// At last, the magnitude.
	sb = math.Sin(b)
	mag = -8.68 + 2.52*math.Abs(up+omg-u) - 2.60*math.Abs(sb) + 1.25*(sb*sb)
	helio()
	geo()
}

func plut() {
	var pturbl, pturbb, pturbr, lograd, dele, enom, vnom, nd, sl,
		capj, capn, eye, comg, omg, sb, su, cu, u, b, up, sd, ca, sa, cy float64
	cy = (eday - elemPlut[0]) / 36525. // Per julian century
	mrad = elemPlut[1] + elemPlut[1+6]*cy
	ecc = elemPlut[2] + elemPlut[2+6]*cy
	cy = cy / 3600 // arcsec/deg per julian century
	incl = elemPlut[3] + elemPlut[3+6]*cy
	node = elemPlut[4] + elemPlut[4+6]*cy
	argp = elemPlut[5] + elemPlut[5+6]*cy
	anom = elemPlut[6] + elemPlut[6+6]*cy - argp
	motion = elemPlut[6+6] / 36525. / 3600
	incl *= radian
	node *= radian
	argp *= radian
	anom = math.Mod(anom, 360.) * radian
	enom = anom + ecc*math.Sin(anom)
	for {
		dele = (anom - enom + ecc*math.Sin(enom)) / (1. - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom = 2. * math.Atan2(math.Sqrt((1.+ecc)/(1.-ecc))*math.Sin(enom/2.), math.Cos(enom/2.))
	rad = mrad * (1. - ecc*math.Cos(enom))
	lambda = vnom + argp
	pturbl = 0.
	lambda += pturbl * radsec
	pturbb = 0.
	pturbr = 0.
	// Reduce to the ecliptic.
	nd = lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl = math.Sin(incl)*math.Sin(nd) + pturbb*radsec
	beta = math.Atan2(sl, pyth(sl))
	lograd = pturbr * 2.30258509
	rad *= 1. + lograd
	lambda -= 1185. * radsec
	beta -= 51. * radsec
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 83.33
	// Computation of magnitude; first, find the geocentric equatorial coordinates of Saturn.
	sd = rad * (math.Cos(beta)*math.Sin(lambda)*math.Sin(obliq) + math.Sin(beta)*math.Cos(obliq))
	sa = rad * (math.Cos(beta)*math.Sin(lambda)*math.Cos(obliq) - math.Sin(beta)*math.Sin(obliq))
	ca = rad * math.Cos(beta) * math.Cos(lambda)
	sd += zms
	sa += yms
	ca += xms
	alpha = math.Atan2(sa, ca)
	delta = math.Atan2(sd, math.Sqrt(sa*sa+ca*ca))
	// Here are the necessary elements of Saturn's rings cf. Exp. Supp. p. 363ff.
	capj = 6.9056 - 0.4322*capt
	capn = 126.3615 + 3.9894*capt + 0.2403*capt2
	eye = 28.0743 - 0.0128*capt
	comg = 168.1179 + 1.3936*capt
	omg = 42.9236 - 2.7390*capt - 0.2344*capt2
	capj *= radian
	capn *= radian
	eye *= radian
	comg *= radian
	omg *= radian
	// Now find saturnicentric ring-plane coords of the earth.
	sb = math.Sin(capj)*math.Cos(delta)*math.Sin(alpha-capn) - math.Cos(capj)*math.Sin(delta)
	su = math.Cos(capj)*math.Cos(delta)*math.Sin(alpha-capn) + math.Sin(capj)*math.Sin(delta)
	cu = math.Cos(delta) * math.Cos(alpha-capn)
	u = math.Atan2(su, cu)
	b = math.Atan2(sb, math.Sqrt(su*su+cu*cu))
	// And then the saturnicentric ring-plane coords of the sun.
	su = math.Sin(eye)*math.Sin(beta) + math.Cos(eye)*math.Cos(beta)*math.Sin(lambda-comg)
	cu = math.Cos(beta) * math.Cos(lambda-comg)
	up = math.Atan2(su, cu)
	// At last, the magnitude.
	sb = math.Sin(b)
	mag = -8.68 + 2.52*math.Abs(up+omg-u) - 2.60*math.Abs(sb) + 1.25*(sb*sb)
	helio()
	geo()
}

type cometElem struct {
	t, // Time of perihelion
	q, // Perihelion distance
	e, // Eccentricity
	i, // Inclination
	w, // Argument of perihelion
	o float64 // Longitude of ascending node
}

func comet() {
	var pturbl, pturbb, pturbr, lograd, dele, enom, vnom, nd, sl float64
	// 153P/Ikeya–Zhang
	t := time.Date(2002, 3, 18, 23, 28, 53, 760000000, time.UTC)
	elem := cometElem{t: timeToJulian(&t) + 2415020, q: 0.5070601, e: 0.990111, i: 28.12106, w: 34.6666, o: 93.1206}
	ecc = elem.e
	if ecc > maxe {
		ecc = maxe
	}
	incl = elem.i * radian
	node = (elem.o + 0.4593) * radian
	argp = (elem.w + elem.o + 0.4066) * radian
	mrad = elem.q / (1 - ecc)
	motion = .01720209895 * math.Sqrt(1/(mrad*mrad*mrad)) / radian
	anom = (eday - (elem.t - 2415020)) * motion * radian
	enom = anom + ecc*math.Sin(anom)
	for {
		dele = (anom - enom + ecc*math.Sin(enom)) / (1. - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom = 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	lambda = vnom + argp
	pturbl = 0
	lambda += pturbl * radsec
	pturbb = 0
	pturbr = 0
	// Reduce to the ecliptic.
	nd = lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl = math.Sin(incl)*math.Sin(nd) + pturbb*radsec
	beta = math.Atan2(sl, math.Sqrt(1-sl*sl))
	lograd = pturbr * 2.30258509
	rad *= 1 + lograd
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 0
	mag = 5.47 + 6.1/2.303*math.Log(rad)
	helio()
	geo()
}

func star() {
	ra = ostar.point[0].ra
	decl2 = ostar.point[0].decl2
	semi2 = ostar.point[0].semi2
	az = ostar.point[0].az
	el = ostar.point[0].el
	mag = ostar.point[0].mag
}

func pstime(d float64) {
	setime(d)
	semi = 0
	motion = 0
	rad = 1e9
	lambda = 0
	beta = 0
	helio()
	geo()
	fmt.Printf(" %f %f %f %4.0f", lha, nlat, awlong, elev/3.28084)
}

func setime(d float64) {
	var x, xm, ym, zm float64
	eday = d + ΔT/86400
	wlong = awlong + 15.*ΔT*radsec
	capt = eday / 36524.220e0
	capt2 = capt * capt
	capt3 = capt * capt2
	nutate()
	eday += .1
	sun()
	srad = rad
	xm = rad * math.Cos(beta) * math.Cos(lambda)
	ym = rad * math.Cos(beta) * math.Sin(lambda)
	zm = rad * math.Sin(beta)
	eday -= .1
	sun()
	xms = rad * math.Cos(beta) * math.Cos(lambda)
	yms = rad * math.Cos(beta) * math.Sin(lambda)
	zms = rad * math.Sin(beta)
	x = .057756
	xdot = x * (xm - xms)
	ydot = x * (ym - yms)
	zdot = x * (zm - zms)
}

func nutate() {
	// Nutation of the equinoxes is a wobble of the pole of the earth's
	// rotation whose magnitude is about 9 seconds of arc and whose period
	// is about 18.6 years. It depends upon the pull of the sun and moon on
	// the equatorial bulge of the earth. phi and eps are the two angles
	// which specify the true pole with respect to the mean pole. All
	// coefficients are from Exp. Supp. pp.44-45.
	var msun, mnom, noded, dmoon float64
	mnom = 296.104608 + 13.0649924465*eday + 9.192e-3*capt2 + 14.38e-6*capt3
	mnom *= radian
	msun = 358.475833 + .9856002669*eday - .150e-3*capt2 - 3.33e-6*capt3
	msun *= radian
	noded = 11.250889 + 13.2293504490*eday - 3.211e-3*capt2 - 0.33e-6*capt3
	noded *= radian
	dmoon = 350.737486 + 12.1907491914*eday - 1.436e-3*capt2 + 1.89e-6*capt3
	dmoon *= radian
	node = 259.183275 - .0529539222*eday + 2.078e-3*capt2 + 2.22e-6*capt3
	node *= radian
	phi = 0.
	eps = 0.
	dphi = 0.
	deps = 0.
	phi = -(17.2327 + .01737*capt) * math.Sin(node)
	phi += sinadd(nutf[0], nutc[0], []float64{node, noded, dmoon, msun})
	eps = cosadd(nutf[1], nutc[1], []float64{node, noded, dmoon, msun})
	dphi = sinadd(nutf[2], nutc[2], []float64{node, noded, mnom, dmoon})
	deps = cosadd(nutf[3], nutc[3], []float64{node, noded, mnom})
	phi = (phi + dphi) * radsec
	eps = (eps + deps) * radsec
	dphi *= radsec
	deps *= radsec
	obliq = 23.452294 - .0130125*capt - 1.64e-6*capt2 + 0.503e-6*capt3
	obliq *= radian
	tobliq = obliq + eps
	gst = 99.690983 + 360.9856473354*eday + .000387*capt2
	gst -= 180.
	gst = math.Mod(gst, 360.)
	if gst < 0. {
		gst += 360.
	}
	gst *= radian
	gst += phi * math.Cos(obliq)
}

func setobj(o *obj1) {
	*o = obj1{ra: ra, decl2: decl2, semi2: semi2, az: az, el: el, mag: mag}
}

func output(n string, p obj1) {
	if n == "" {
		fmt.Printf(" SAO %s", sao)
	} else {
		fmt.Printf("%10s", n)
	}
	fmt.Printf(" %f %f %9.4f %9.4f %9.4f", p.ra, p.decl2, p.az, p.el, p.semi2)
	if n == osun.name || n == omoon.name {
		fmt.Printf(" %7.4f", p.mag)
	}
	fmt.Println()
}

func dist(o1, o2 obj1) float64 {
	d := math.Sin(o1.decl2)*math.Sin(o2.decl2) +
		math.Cos(o1.decl2)*math.Cos(o2.decl2)*math.Cos(o1.ra-o2.ra)
	return math.Abs(math.Atan2(pyth(d), d)) / radsec
}

func search() error {
	for i, o := range objs {
		if o.name == oshad.name {
			continue
		}
		t := rise(*o, -.833)
		var flag int
		if i == 0 {
			flag = ptime
		} else {
			flag = ptime | dark
		}
		if t >= 0 {
			err := event(evt{s: fmt.Sprintf("%s rises at ", o.fname), tim: t, flag: flag})
			if err != nil {
				return err
			}
		}
		t = set(*o, -.833)
		if t >= 0 {
			err := event(evt{s: fmt.Sprintf("%s sets at ", o.fname), tim: t, flag: flag})
			if err != nil {
				return err
			}
		}
		if o.name == osun.name {
			for j := range 4 {
				t = solstice(j)
				if t >= 0 {
					err := event(evt{s: fmt.Sprintf("%s at ", solstr[j]), tim: t, flag: signif | ptime})
					if err != nil {
						return err
					}
				}
			}
			bettab := []struct {
				beta   float64
				shower string
			}{
				{-1.3572, "Quadrantid"},
				{0.7620, "Eta aquarid"},
				{1.5497, "Ophiuchid"},
				{2.1324, "Capricornid"},
				{2.1991, "Delta aquarid"},
				{2.2158, "Pisces australid"},
				{2.4331, "Perseid"},
				{-2.6578, "Orionid"},
				{-1.8678, "Phoenicid"},
				{-1.7260, "Geminid"},
			}
			for _, b := range bettab {
				t = betcross(b.beta)
				if t >= 0 {
					err := event(evt{s: fmt.Sprintf("%s meteor shower", b.shower), tim: t, flag: signif})
					if err != nil {
						return err
					}
				}
			}
			t = rise(*o, -18)
			if t >= 0 {
				err := event(evt{s: "Twilight starts at ", tim: t, flag: ptime})
				if err != nil {
					return err
				}
			}
			t = set(*o, -18)
			if t >= 0 {
				err := event(evt{s: "Twilight ends at ", tim: t, flag: ptime})
				if err != nil {
					return err
				}
			}
		}
		if o.name == omoon.name {
			for j := range len(o.point) - 2 {
				if o.point[j].mag > .75 && o.point[j+1].mag < .25 {
					err := event(evt{s: "New moon"})
					if err != nil {
						return err
					}
				}
				if o.point[j].mag <= .25 && o.point[j+1].mag > .25 {
					err := event(evt{s: "First quarter moon"})
					if err != nil {
						return err
					}
				}
				if o.point[j].mag <= .5 && o.point[j+1].mag > .5 {
					err := event(evt{s: "Full moon"})
					if err != nil {
						return err
					}
				}
				if o.point[j].mag <= .75 && o.point[j+1].mag > .75 {
					err := event(evt{s: "Last quarter moon"})
					if err != nil {
						return err
					}
				}
			}
		}
		if o.name == omerc.name || o.name == ovenus.name {
			t = float64(melong(*o))
			if t >= 0 {
				t = rise(*o, 0) - rise(osun, 0)
				if t < 0 {
					t += npts
				}
				if t > npts {
					t -= npts
				}
				if t > npts/2 {
					err := event(evt{s: fmt.Sprintf("Morning elongation of %s", o.fname), flag: signif})
					if err != nil {
						return err
					}
				} else {
					err := event(evt{s: fmt.Sprintf("Evening elongation of %s", o.fname), flag: signif})
					if err != nil {
						return err
					}
				}
			}
		}
		for _, p := range objs[i+1:] {
			if o.name == omoon.name || p.name == omoon.name {
				if err := occult(*o, *p); err != nil {
					return err
				}
				if occ.t3 < 0 {
					continue
				}
				if o.name == osun.name || p.name == omoon.name {
					if occ.t1 >= 0 {
						err := event(evt{s: fmt.Sprintf("Partial eclipse of %s begins at ", o.fname), tim: occ.t1, flag: signif | ptime})
						if err != nil {
							return err
						}
					}
					if occ.t2 >= 0 {
						err := event(evt{s: fmt.Sprintf("Total eclipse of %s begins at ", o.fname), tim: occ.t2, flag: signif | ptime})
						if err != nil {
							return err
						}
					}
					if occ.t4 >= 0 {
						err := event(evt{s: fmt.Sprintf("Total eclipse of %s ends at ", o.fname), tim: occ.t4, flag: signif | ptime})
						if err != nil {
							return err
						}
					}
					if occ.t5 >= 0 {
						err := event(evt{s: fmt.Sprintf("Partial eclipse of %s ends at ", o.fname), tim: occ.t5, flag: signif | ptime})
						if err != nil {
							return err
						}
					}
				} else {
					if occ.t1 >= 0 {
						err := event(evt{s: fmt.Sprintf("Occultation of %s begins at ", o.fname), tim: occ.t1, flag: signif | ptime})
						if err != nil {
							return err
						}
					}
					if occ.t5 >= 0 {
						err := event(evt{s: fmt.Sprintf("Occultation of %s ends at ", o.fname), tim: occ.t5, flag: signif | ptime})
						if err != nil {
							return err
						}
					}
				}
				continue
			}
			if o.name == osun.name {
				if p.name != omerc.name && p.name != ovenus.name {
					continue
				}
				if err := occult(*o, *p); err != nil {
					return err
				}
				if occ.t3 >= 0 {
					if occ.t1 >= 0 {
						err := event(evt{s: fmt.Sprintf("Transit of %s begins at ", p.fname), tim: occ.t1, flag: signif | light | ptime})
						if err != nil {
							return err
						}
					}
					if occ.t5 >= 0 {
						err := event(evt{s: fmt.Sprintf("Transit of %s ends at ", p.fname), tim: occ.t5, flag: signif | light | ptime})
						if err != nil {
							return err
						}
					}
				}
				continue
			}
			t = dist(o.point[0], p.point[0])
			if t > 5000 {
				continue
			}
			err := event(evt{s: fmt.Sprintf("%s is in the house of %s", o.fname, p.fname)})
			if err != nil {
				return err
			}
		}
	}
	if *oflag {
		if err := stars(); err != nil {
			return err
		}
	}
	evflush()
	return nil
}

func rise(o obj2, el float64) float64 {
	for i := 1; i < len(o.point); i++ {
		e1, e2 := o.point[i-1].el, o.point[i].el
		if e1 <= el && e2 > el {
			return float64(i-1) + (el-e1)/(e2-e1)
		}
	}
	return -1
}

func set(o obj2, el float64) float64 {
	for i := 1; i < len(o.point); i++ {
		e1, e2 := o.point[i-1].el, o.point[i].el
		if e1 > el && e2 <= el {
			return float64(i-1) + (el-e1)/(e2-e1)
		}
	}
	return -1
}

func solstice(n int) float64 {
	d3 := float64(n)*math.Pi/2 - math.Pi
	if n == 0 {
		d3 += math.Pi
	}
	var d1, d2 float64
	for i := range len(osun.point) - 1 {
		d1, d2 = d2, osun.point[i].ra
		if n == 0 {
			d2 -= math.Pi
			if d2 < -math.Pi {
				d2 += pipi
			}
		}
		if i >= 1 && d3 >= d1 && d3 < d2 {
			return float64(i) - (d3-d2)/(d1-d2)
		}
	}
	return -1
}

func betcross(b float64) float64 {
	for i := 1; i < len(osun.point); i++ {
		d1, d2 := osun.point[i-1].mag, osun.point[i].mag
		if b >= d1 && b < d2 {
			return float64(i) - (b-d2)/(d1-d2)
		}
	}
	return -1
}

func melong(o obj2) int {
	for i := 2; i < len(o.point); i++ {
		d1 := dist(o.point[i-2], osun.point[i-2])
		d2 := dist(o.point[i-1], osun.point[i-1])
		d3 := dist(o.point[i], osun.point[i])
		if d2 >= d1 && d2 >= d3 {
			return i - 2
		}
	}
	return -1
}

func event(e evt) error {
	if e.flag&dark > 0 && sunel(e.tim) > -12 {
		return nil
	}
	if e.flag&light > 0 && sunel(e.tim) < 0 {
		return nil
	}
	if len(events) >= nevents {
		return errors.New("too many events")
	}
	events = append(events, e)
	return nil
}

func evflush() {
	slices.SortFunc(events, func(e1, e2 evt) int {
		t1, t2 := e1.tim, e2.tim
		if e1.flag&signif > 0 {
			t1 -= 1000.
		}
		if e2.flag&signif > 0 {
			t2 -= 1000.
		}
		return cmp.Compare(t1, t2)
	})
	for _, e := range events {
		if e.flag&ptime > 0 {
			fmt.Printf("%s%s\n", e.s, julianToTime(day+e.tim*deld))
		} else {
			fmt.Printf("%s\n", e.s)
		}
	}
}

func sunel(t float64) float64 {
	i := int(t)
	if i < 0 || i > npts {
		return -90
	}
	return osun.point[i].el + (t-float64(i))*(osun.point[i+1].el-osun.point[i].el)
}

func occult(o1, o2 obj2) error {
	occ.t1 = -100
	occ.t2 = -100
	occ.t3 = -100
	occ.t4 = -100
	occ.t5 = -100
	var i int
	var d1, d2, d3 float64
	var ok bool
	for i = 2; i < len(o1.point); i++ {
		d1 = dist(o1.point[i-2], o2.point[i-2])
		d2 = dist(o1.point[i-1], o2.point[i-1])
		d3 = dist(o1.point[i], o2.point[i])
		if d2 <= d1 && d2 <= d3 {
			ok = true
			break
		}
	}
	if !ok {
		return nil
	}
	ok = false
	n := 2880 * per / npts // 1 min steps
	i -= 2
	set3pt(o1, i, &occ1)
	set3pt(o2, i, &occ2)
	di := float64(i)
	x := 0.
	dx := 2 / n
	for i = range int(n + 1) {
		setpt(&occ1, x)
		setpt(&occ2, x)
		d1, d2 = d2, d3
		d3 = dist(occ1.act, occ2.act)
		if i >= 2 && d2 <= d1 && d2 <= d3 {
			ok = true
			break
		}
		x += dx
	}
	if !ok {
		return errors.New("bad 1 \n")
	}
	ok = false
	if d2 > occ1.act.semi2+occ2.act.semi2+50 {
		return nil
	}
	di += x - 3*dx
	x = 0
	var xo1, xo2 obj2
	for i = range 3 {
		setime(day + deld*(di+x))
		o1.f()
		setobj(&xo1.point[i])
		o2.f()
		setobj(&xo2.point[i])
		x += 2 * dx
	}
	dx /= 60
	x = 0
	set3pt(xo1, 0, &occ1)
	set3pt(xo2, 0, &occ2)
	for i = range 241 {
		setpt(&occ1, x)
		setpt(&occ2, x)
		d1, d2 = d2, d3
		d3 = dist(occ1.act, occ2.act)
		if i >= 2 && d2 <= d1 && d2 <= d3 {
			ok = true
			break
		}
		x += 1. / 120
	}
	if !ok {
		return errors.New("bad 2 \n")
	}
	if d2 > occ1.act.semi2+occ2.act.semi2 {
		return nil
	}
	i1 := i - 1
	x1 := x - 1./120
	occ.t3 = di + float64(i1)*dx
	occ.e3 = occ1.act.el
	d3 = occ1.act.semi2 - occ2.act.semi2
	if d3 < 0 {
		d3 = -d3
	}
	d4 := occ1.act.semi2 + occ2.act.semi2
	i = i1
	x = x1
	for {
		setpt(&occ1, x)
		setpt(&occ2, x)
		d1 = d2
		d2 = dist(occ1.act, occ2.act)
		if i != i1 {
			if d1 <= d3 && d2 > d3 {
				occ.t4 = di + (float64(i)-.5)*dx
				occ.e4 = occ1.act.el
			}
			if d2 > d4 {
				if d1 <= d4 {
					occ.t5 = di + (float64(i)-.5)*dx
					occ.e5 = occ1.act.el
				}
				break
			}
		}
		x += 1. / 120
		i++
	}
	i = i1
	x = x1
	for {
		setpt(&occ1, x)
		setpt(&occ2, x)
		d1 = d2
		d2 = dist(occ1.act, occ2.act)
		if i != i1 {
			if d1 <= d3 && d2 > d3 {
				occ.t2 = di + (float64(i)-.5)*dx
				occ.e2 = occ1.act.el
			}
			if d2 > d4 {
				if d1 <= d4 {
					occ.t1 = di + (float64(i)-.5)*dx
					occ.e1 = occ1.act.el
				}
				break
			}
		}
		x -= 1. / 120
		i--
	}
	return nil
}

func set3pt(o obj2, i int, oc *occt) {
	p1, p2, p3 := o.point[i], o.point[i+1], o.point[i+2]
	oc.del0.ra = p1.ra
	oc.del0.decl2 = p1.decl2
	oc.del0.semi2 = p1.semi2
	oc.del0.el = p1.el
	a := p2.ra - p1.ra
	oc.del1.ra = pinorm(a)
	a = p2.decl2 - p1.decl2
	oc.del1.decl2 = pinorm(a)
	oc.del1.semi2 = p2.semi2 - p1.semi2
	oc.del1.el = p2.el - p1.el
	a = p1.ra + p3.ra - 2*p2.ra
	oc.del2.ra = pinorm(a) / 2
	a = p1.decl2 + p3.decl2 - 2*p2.decl2
	oc.del2.decl2 = pinorm(a) / 2
	oc.del2.semi2 = (p1.semi2 + p3.semi2 - 2*p2.semi2) / 2
	oc.del2.el = (p1.el + p3.el - 2*p2.el) / 2
}

func setpt(o *occt, x float64) {
	y := x * (x - 1)
	o.act.ra = o.del0.ra + x*o.del1.ra + y*o.del2.ra
	o.act.decl2 = o.del0.decl2 + x*o.del1.decl2 + y*o.del2.decl2
	o.act.semi2 = o.del0.semi2 + x*o.del1.semi2 + y*o.del2.semi2
	o.act.el = o.del0.el + x*o.del1.el + y*o.del2.el
}

func stars() error {
	sd := 1000 * radsec
	lomoon := omoon.point[0].ra - sd
	if lomoon < 0 {
		lomoon += pipi
	}
	himoon := omoon.point[npts+1].ra + sd
	if himoon > pipi {
		himoon -= pipi
	}
	lomoon *= 12 / math.Pi
	himoon *= 12 / math.Pi
	wrap := false
	if lomoon > himoon {
		wrap = true
	}
	f, err := os.Open(filepath.Join(root, "sky", "estartab"))
	if err != nil {
		return err
	}
	defer f.Close()
	epoch := (1950-1900)*365.24220 + 0.313
	s := bufio.NewScanner(f)
	for s.Scan() {
		l := s.Text()
		// Read mean places of stars at epoch of star table.
		rah, err := strconv.Atoi(l[18:20])
		if err != nil {
			return err
		}
		ram, err := strconv.Atoi(l[21:23])
		if err != nil {
			return err
		}
		ras, err := strconv.ParseFloat(l[24:30], 64)
		if err != nil {
			return err
		}
		alpha = float64(rah) + float64(ram)/60 + ras/3600
		if wrap && alpha < lomoon && alpha > himoon {
			continue
		} else if alpha < lomoon || alpha > himoon {
			continue
		}
		sao = l[:6]
		da, err := strconv.ParseFloat(strings.TrimSpace(l[31:37]), 64)
		if err != nil {
			return err
		}
		dday, err := strconv.Atoi(strings.TrimSpace(l[38:41]))
		if err != nil {
			return err
		}
		dmin, err := strconv.Atoi(l[42:44])
		if err != nil {
			return err
		}
		dsec, err := strconv.ParseFloat(l[45:50], 64)
		if err != nil {
			return err
		}
		dd, err := strconv.ParseFloat(strings.TrimSpace(l[51:57]), 64)
		if err != nil {
			return err
		}
		px, err := strconv.ParseFloat(l[58:61], 64)
		if err != nil {
			return err
		}
		mag, err := strconv.ParseFloat(l[63:67], 64)
		if err != nil {
			return err
		}
		// Convert rt ascension and declination to internal format.
		delta = float64(abs(dday)) + float64(dmin)/60 + dsec/3600
		if dday < 0 {
			delta = -delta
		}
		// Remove E-terms of aberration except when finding catalog mean places.
		alpha += (.341 / (3600. * 15.)) * math.Sin((alpha+11.26)*15.*radian) / math.Cos(delta*radian)
		delta += (.341/3600.)*math.Cos((alpha+11.26)*15.*radian)*math.Sin(delta*radian) - (.029/3600.)*math.Cos(delta*radian)
		// Correct for proper motion.
		tau := (eday - epoch) / 365.24220
		alpha += tau * da / 3600.
		delta += tau * dd / 3600.
		alpha *= 15. * radian
		delta *= radian
		// Convert to rectangular coordinates merely for convenience.
		xm := math.Cos(delta) * math.Cos(alpha)
		ym := math.Cos(delta) * math.Sin(alpha)
		zm := math.Sin(delta)
		// Convert mean places at epoch of startable to current epoch (i.e. compute relevant precession).
		capt0 := (epoch - 18262.427) / 36524.220e0
		capt1 := (eday - epoch) / 36524.220
		capt12 := capt1 * capt1
		capt13 := capt12 * capt1
		xx := -(.00029696+26.e-8*capt0)*capt12 - 13.e-8*capt13
		yx := -(.02234941+1355.e-8*capt0)*capt1 - 676.e-8*capt12 + 221.e-8*capt13
		zx := -(.00971690-414.e-8*capt0)*capt1 + 207.e-8*capt12 + 96.e-8*capt13
		yy := -(.00024975+30.e-8*capt0)*capt12 - 15.e-8*capt13
		zy := -(.00010858 + 2.e-8*capt0) * capt12
		zz := -(.00004721 - 4.e-8*capt0) * capt12
		dxm := xx*xm + yx*ym + zx*zm
		dym := -yx*xm + yy*ym + zy*zm
		dzm := -zx*xm + zy*ym + zz*zm
		xm += dxm
		ym += dym
		zm += dzm
		// Convert to mean ecliptic system of date.
		alpha = math.Atan2(ym, xm)
		delta = math.Atan2(zm, math.Sqrt(xm*xm+ym*ym))
		cl := math.Cos(delta) * math.Cos(alpha)
		sl := math.Cos(delta)*math.Sin(alpha)*math.Cos(obliq) + math.Sin(delta)*math.Sin(obliq)
		sb := -math.Cos(delta)*math.Sin(alpha)*math.Sin(obliq) + math.Sin(delta)*math.Cos(obliq)
		lambda = math.Atan2(sl, cl)
		beta = math.Atan2(sb, math.Sqrt(cl*cl+sl*sl))
		rad = 1.e9
		if px != 0 {
			rad = 20600 / px
		}
		motion = 0
		semi = 0
		helio()
		geo()
		sd = .0896833e0*math.Cos(beta)*math.Sin(lambda-1.3820+.00092422117*eday) + 0.99597*math.Sin(beta)
		if math.Abs(sd) > .0183 {
			continue
		}
		for _, p := range ostar.point {
			setobj(&p)
		}
		if err = occult(omoon, ostar); err != nil {
			return err
		}
		if occ.t1 >= 0 || occ.t5 >= 0 {
			i := ptime
			if mag > 2 {
				i |= dark
			}
			if mag < 5 {
				i |= signif
			}
			if occ.t1 >= 0 && occ.e1 >= 0 {
				err := event(evt{s: fmt.Sprintf("Occultation of SAO %s begins at ", sao), tim: occ.t1, flag: i})
				if err != nil {
					return err
				}
			}
			if occ.t5 >= 0 && occ.e5 >= 0 {
				err := event(evt{s: fmt.Sprintf("Occultation of SAO %s ends at ", sao), tim: occ.t5, flag: i})
				if err != nil {
					return err
				}
			}
		}
	}
	return nil
}
