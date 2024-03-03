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
	secondsPerDay = 24 * 60 * 60
	twoPi         = 2 * math.Pi
	radian        = math.Pi / 180
	radsec        = radian / 3600
	converge      = 1e-14
	metersToFeet  = 3.28084
	iVal          = 1.0
	numPoints     = 12
	stepSize      = iVal / numPoints
	dark          = 1 << iota
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
	point       [numPoints + 2]obj1
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
	julian       = flag.Bool("j", false, "print Julian date")
	printPos     = flag.Bool("p", false, "print positions of objects at the given time")
	searchOccult = flag.Bool("o", false, "search for stellar occultations")
	local        = flag.Bool("k", false, "print times in local time")
	includeComet = flag.Bool("m", false, "include single comet in the list of objects")
	periods      = flag.Int("c", 1, "report for n successive days")
	interval     = flag.Float64("C", iVal, "used with -c, set the interval to d days")
	startDate    = flag.String("d", "", "read start date")
	eclipse      = flag.String("e", "", "report distance between the centers of objects")
	loc          = flag.String("l", "", "read latitude, longitude, and elevation")
	dt           = flag.Float64("t", 0, "read ΔT")

	root = os.Getenv("PLAN9")
	wlong, awlong, nlat, elev,
	obliq, phi, eps, tobliq,
	day, eday, capt, capt2, capt3, gst,
	ΔT, erad, glat,
	xms, yms, zms, xdot, ydot, zdot,
	motion, lambda, beta, rad, mag, semi,
	alpha, delta, hp,
	ra, semi2, lha, decl2, lmb2, az, el,
	meday, seday, mhp, salph, sdelt, srad float64
	sao    string
	events []evt
	oSun   = obj2{name: "sun", fname: "The sun", f: fSun}
	oMoon  = obj2{name: "moon", fname: "The moon", f: moon}
	oShad  = obj2{name: "shadow", fname: "The shadow", f: shad}
	oMerc  = obj2{name: "mercury", fname: "Mercury", f: merc}
	oVenus = obj2{name: "venus", fname: "Venus", f: venus}
	objs   = []*obj2{
		&oSun, &oMoon, &oShad, &oMerc, &oVenus,
		{name: "mars", fname: "Mars", f: mars},
		{name: "jupiter", fname: "Jupiter", f: jup},
		{name: "saturn", fname: "Saturn", f: sat},
		{name: "uranus", fname: "Uranus", f: uran},
		{name: "neptune", fname: "Neptune", f: nept},
		{name: "pluto", fname: "Pluto", f: plut},
		{name: "comet", fname: "Comet", f: comet},
	}
	oStar      obj2
	occ        obj3
	occ1, occ2 occt
	moonTabs   = [...][]moonTab{
		{
			{f: 0.127, c: [4]int{0, 0, 0, 6}},
			{f: 13.902, c: [4]int{0, 0, 0, 4}},
			{f: 2369.912, c: [4]int{0, 0, 0, 2}},
			{f: 1.979, c: [4]int{1, 0, 0, 4}},
			{f: 191.953, c: [4]int{1, 0, 0, 2}},
			{f: 22639.5, c: [4]int{1, 0, 0, 0}},
			{f: -4586.465, c: [4]int{1, 0, 0, -2}},
			{f: -38.428, c: [4]int{1, 0, 0, -4}},
			{f: -0.393, c: [4]int{1, 0, 0, -6}},
			{f: -0.289, c: [4]int{0, 1, 0, 4}},
			{f: -24.42, c: [4]int{0, 1, 0, 2}},
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
			{f: -0.57, c: [4]int{2, 0, 0, -6}},
			{f: -2.921, c: [4]int{1, 1, 0, 2}},
			{f: -109.673, c: [4]int{1, 1, 0, 0}},
			{f: -205.962, c: [4]int{1, 1, 0, -2}},
			{f: -4.391, c: [4]int{1, 1, 0, -4}},
			{f: -0.072, c: [4]int{1, 1, 0, -6}},
			{f: 0.283, c: [4]int{1, -1, 0, 4}},
			{f: 14.577, c: [4]int{1, -1, 0, 2}},
			{f: 147.687, c: [4]int{1, -1, 0, 0}},
			{f: 28.475, c: [4]int{1, -1, 0, -2}},
			{f: 0.636, c: [4]int{1, -1, 0, -4}},
			{f: -0.189, c: [4]int{0, 2, 0, 2}},
			{f: -7.486, c: [4]int{0, 2, 0, 0}},
			{f: -8.096, c: [4]int{0, 2, 0, -2}},
			{f: -0.151, c: [4]int{0, 2, 0, -4}},
			{f: -0.085, c: [4]int{0, 0, 2, 4}},
			{f: -5.741, c: [4]int{0, 0, 2, 2}},
			{f: -411.608, c: [4]int{0, 0, 2, 0}},
			{f: -55.173, c: [4]int{0, 0, 2, -2}},
			{f: -8.466, c: [4]int{1, 0, 0, 1}},
			{f: 18.609, c: [4]int{1, 0, 0, -1}},
			{f: 3.215, c: [4]int{1, 0, 0, -3}},
			{f: 0.15, c: [4]int{0, 1, 0, 3}},
			{f: 18.023, c: [4]int{0, 1, 0, 1}},
			{f: 0.56, c: [4]int{0, 1, 0, -1}},
			{f: 1.06, c: [4]int{3, 0, 0, 2}},
			{f: 36.124, c: [4]int{3, 0, 0, 0}},
			{f: -13.193, c: [4]int{3, 0, 0, -2}},
			{f: -1.187, c: [4]int{3, 0, 0, -4}},
			{f: -0.293, c: [4]int{3, 0, 0, -6}},
			{f: -0.29, c: [4]int{2, 1, 0, 2}},
			{f: -7.649, c: [4]int{2, 1, 0, 0}},
			{f: -8.627, c: [4]int{2, 1, 0, -2}},
			{f: -2.74, c: [4]int{2, 1, 0, -4}},
			{f: -0.091, c: [4]int{2, 1, 0, -6}},
			{f: 1.181, c: [4]int{2, -1, 0, 2}},
			{f: 9.703, c: [4]int{2, -1, 0, 0}},
			{f: -2.494, c: [4]int{2, -1, 0, -2}},
			{f: 0.36, c: [4]int{2, -1, 0, -4}},
			{f: -1.167, c: [4]int{1, 2, 0, 0}},
			{f: -7.412, c: [4]int{1, 2, 0, -2}},
			{f: -0.311, c: [4]int{1, 2, 0, -4}},
			{f: 0.757, c: [4]int{1, -2, 0, 2}},
			{f: 2.58, c: [4]int{1, -2, 0, 0}},
			{f: 2.533, c: [4]int{1, -2, 0, -2}},
			{f: -0.103, c: [4]int{0, 3, 0, 0}},
			{f: -0.344, c: [4]int{0, 3, 0, -2}},
			{f: -0.992, c: [4]int{1, 0, 2, 2}},
			{f: -45.099, c: [4]int{1, 0, 2, 0}},
			{f: -0.179, c: [4]int{1, 0, 2, -2}},
			{f: -0.301, c: [4]int{1, 0, 2, -4}},
			{f: -6.382, c: [4]int{1, 0, -2, 2}},
			{f: 39.528, c: [4]int{1, 0, -2, 0}},
			{f: 9.366, c: [4]int{1, 0, -2, -2}},
			{f: 0.202, c: [4]int{1, 0, -2, -4}},
			{f: 0.415, c: [4]int{0, 1, 2, 0}},
			{f: -2.152, c: [4]int{0, 1, 2, -2}},
			{f: -1.44, c: [4]int{0, 1, -2, 2}},
			{f: 0.076, c: [4]int{0, 1, -2, 0}},
			{f: 0.384, c: [4]int{0, 1, -2, -2}},
			{f: -0.586, c: [4]int{2, 0, 0, 1}},
			{f: 1.75, c: [4]int{2, 0, 0, -1}},
			{f: 1.225, c: [4]int{2, 0, 0, -3}},
			{f: 1.267, c: [4]int{1, 1, 0, 1}},
			{f: 0.137, c: [4]int{1, 1, 0, -1}},
			{f: 0.233, c: [4]int{1, 1, 0, -3}},
			{f: -0.122, c: [4]int{1, -1, 0, 1}},
			{f: -1.089, c: [4]int{1, -1, 0, -1}},
			{f: -0.276, c: [4]int{1, -1, 0, -3}},
			{f: 0.255, c: [4]int{0, 0, 2, 1}},
			{f: 0.584, c: [4]int{0, 0, 2, -1}},
			{f: 0.254, c: [4]int{0, 0, 2, -3}},
			{f: 0.07, c: [4]int{4, 0, 0, 2}},
			{f: 1.938, c: [4]int{4, 0, 0, 0}},
			{f: -0.952, c: [4]int{4, 0, 0, -2}},
			{f: -0.551, c: [4]int{3, 1, 0, 0}},
			{f: -0.482, c: [4]int{3, 1, 0, -2}},
			{f: -0.1, c: [4]int{3, 1, 0, -4}},
			{f: 0.088, c: [4]int{3, -1, 0, 2}},
			{f: 0.681, c: [4]int{3, -1, 0, 0}},
			{f: -0.183, c: [4]int{3, -1, 0, -2}},
			{f: -0.297, c: [4]int{2, 2, 0, -2}},
			{f: -0.161, c: [4]int{2, 2, 0, -4}},
			{f: 0.197, c: [4]int{2, -2, 0, 0}},
			{f: 0.254, c: [4]int{2, -2, 0, -2}},
			{f: -0.25, c: [4]int{1, 3, 0, -2}},
			{f: -0.123, c: [4]int{2, 0, 2, 2}},
			{f: -3.996, c: [4]int{2, 0, 2, 0}},
			{f: 0.557, c: [4]int{2, 0, 2, -2}},
			{f: -0.459, c: [4]int{2, 0, -2, 2}},
			{f: -1.37, c: [4]int{2, 0, -2, 0}},
			{f: 0.538, c: [4]int{2, 0, -2, -2}},
			{f: 0.173, c: [4]int{2, 0, -2, -4}},
			{f: 0.263, c: [4]int{1, 1, 2, 0}},
			{f: 0.083, c: [4]int{1, 1, -2, 2}},
			{f: -0.083, c: [4]int{1, 1, -2, 0}},
			{f: 0.426, c: [4]int{1, 1, -2, -2}},
			{f: -0.304, c: [4]int{1, -1, 2, 0}},
			{f: -0.372, c: [4]int{1, -1, -2, 2}},
			{f: 0.083, c: [4]int{1, -1, -2, 0}},
			{f: 0.418, c: [4]int{0, 0, 4, 0}},
			{f: 0.074, c: [4]int{0, 0, 4, -2}},
			{f: 0.13, c: [4]int{3, 0, 0, -1}},
			{f: 0.092, c: [4]int{2, 1, 0, 1}},
			{f: 0.084, c: [4]int{2, 1, 0, -3}},
			{f: -0.352, c: [4]int{2, -1, 0, -1}},
			{f: 0.113, c: [4]int{5, 0, 0, 0}},
			{f: -0.33, c: [4]int{3, 0, 2, 0}},
			{f: 0.09, c: [4]int{1, 0, 4, 0}},
			{f: -0.08, c: [4]int{1, 0, -4, 0}},
		},
		{
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
			{f: -16.4, c: [4]int{3, 0, 0, -2}},
			{f: 3.6, c: [4]int{4, 0, 0, 0}},
			{f: -1.58, c: [4]int{4, 0, 0, -2}},
			{f: -1.59, c: [4]int{0, 1, 0, 4}},
			{f: -25.1, c: [4]int{0, 1, 0, 2}},
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
			{f: -31.7, c: [4]int{-1, 1, 0, -2}},
			{f: -1.53, c: [4]int{-1, 1, 0, -4}},
			{f: -10.56, c: [4]int{2, 1, 0, 0}},
			{f: -7.59, c: [4]int{2, 1, 0, -2}},
			{f: -2.54, c: [4]int{2, 1, 0, -4}},
			{f: 3.32, c: [4]int{2, -1, 0, 2}},
			{f: 11.67, c: [4]int{2, -1, 0, 0}},
			{f: -6.12, c: [4]int{1, 2, 0, -2}},
			{f: -2.4, c: [4]int{-1, 2, 0, 2}},
			{f: -2.32, c: [4]int{-1, 2, 0, 0}},
			{f: -1.82, c: [4]int{-1, 2, 0, -2}},
			{f: -52.14, c: [4]int{0, 0, 2, -2}},
			{f: -1.67, c: [4]int{0, 0, 2, -4}},
			{f: -9.52, c: [4]int{1, 0, 2, -2}},
			{f: -85.13, c: [4]int{-1, 0, 2, 0}},
			{f: 3.37, c: [4]int{-1, 0, 2, -2}},
			{f: -2.26, c: [4]int{0, 1, 2, -2}},
		},
		{
			{f: -0.725, c: [4]int{0, 0, 0, 1}},
			{f: 0.601, c: [4]int{0, 0, 0, 2}},
			{f: 0.394, c: [4]int{0, 0, 0, 3}},
			{f: -0.445, c: [4]int{1, 0, 0, 4}},
			{f: 0.455, c: [4]int{1, 0, 0, 1}},
			{f: 0.192, c: [4]int{1, 0, 0, -3}},
			{f: 5.679, c: [4]int{2, 0, 0, -2}},
			{f: -0.308, c: [4]int{2, 0, 0, -4}},
			{f: -0.166, c: [4]int{3, 0, 0, 2}},
			{f: -1.3, c: [4]int{3, 0, 0, 0}},
			{f: 0.258, c: [4]int{3, 0, 0, -2}},
			{f: -1.302, c: [4]int{0, 1, 0, 0}},
			{f: -0.416, c: [4]int{0, 1, 0, -4}},
			{f: -0.74, c: [4]int{0, 2, 0, -2}},
			{f: 0.787, c: [4]int{1, 1, 0, 2}},
			{f: 0.461, c: [4]int{1, 1, 0, 0}},
			{f: 2.056, c: [4]int{1, 1, 0, -2}},
			{f: -0.471, c: [4]int{1, 1, 0, -4}},
			{f: -0.443, c: [4]int{-1, 1, 0, 2}},
			{f: 0.679, c: [4]int{-1, 1, 0, 0}},
			{f: -1.54, c: [4]int{-1, 1, 0, -2}},
			{f: 0.259, c: [4]int{2, 1, 0, 0}},
			{f: -0.212, c: [4]int{2, -1, 0, 2}},
			{f: -0.151, c: [4]int{2, -1, 0, 0}},
		},
		{
			{f: -526.069, c: [4]int{0, 0, 1, -2}},
			{f: -3.352, c: [4]int{0, 0, 1, -4}},
			{f: 44.297, c: [4]int{1, 0, 1, -2}},
			{f: -6, c: [4]int{1, 0, 1, -4}},
			{f: 20.599, c: [4]int{-1, 0, 1, 0}},
			{f: -30.598, c: [4]int{-1, 0, 1, -2}},
			{f: -24.649, c: [4]int{-2, 0, 1, 0}},
			{f: -2, c: [4]int{-2, 0, 1, -2}},
			{f: -22.571, c: [4]int{0, 1, 1, -2}},
			{f: 10.985, c: [4]int{0, -1, 1, -2}},
		},
		{
			{f: 0.2607, c: [4]int{0, 0, 0, 4}},
			{f: 28.2333, c: [4]int{0, 0, 0, 2}},
			{f: 0.0433, c: [4]int{1, 0, 0, 4}},
			{f: 3.0861, c: [4]int{1, 0, 0, 2}},
			{f: 186.5398, c: [4]int{1, 0, 0, 0}},
			{f: 34.3117, c: [4]int{1, 0, 0, -2}},
			{f: 0.6008, c: [4]int{1, 0, 0, -4}},
			{f: -0.3, c: [4]int{0, 1, 0, 2}},
			{f: -0.3997, c: [4]int{0, 1, 0, 0}},
			{f: 1.9178, c: [4]int{0, 1, 0, -2}},
			{f: 0.0339, c: [4]int{0, 1, 0, -4}},
			{f: -0.9781, c: [4]int{0, 0, 0, 1}},
			{f: 0.2833, c: [4]int{2, 0, 0, 2}},
			{f: 10.1657, c: [4]int{2, 0, 0, 0}},
			{f: -0.3039, c: [4]int{2, 0, 0, -2}},
			{f: 0.3722, c: [4]int{2, 0, 0, -4}},
			{f: 0.0109, c: [4]int{2, 0, 0, -6}},
			{f: -0.0484, c: [4]int{1, 1, 0, 2}},
			{f: -0.949, c: [4]int{1, 1, 0, 0}},
			{f: 1.4437, c: [4]int{1, 1, 0, -2}},
			{f: 0.0673, c: [4]int{1, 1, 0, -4}},
			{f: 0.2302, c: [4]int{1, -1, 0, 2}},
			{f: 1.1528, c: [4]int{1, -1, 0, 0}},
			{f: -0.2257, c: [4]int{1, -1, 0, -2}},
			{f: -0.0102, c: [4]int{1, -1, 0, -4}},
			{f: 0.0918, c: [4]int{0, 2, 0, -2}},
			{f: -0.0124, c: [4]int{0, 0, 2, 0}},
			{f: -0.1052, c: [4]int{0, 0, 2, -2}},
			{f: -0.1093, c: [4]int{1, 0, 0, 1}},
			{f: 0.0118, c: [4]int{1, 0, 0, -1}},
			{f: -0.0386, c: [4]int{1, 0, 0, -3}},
			{f: 0.1494, c: [4]int{0, 1, 0, 1}},
			{f: 0.0243, c: [4]int{3, 0, 0, 2}},
			{f: 0.6215, c: [4]int{3, 0, 0, 0}},
			{f: -0.1187, c: [4]int{3, 0, 0, -2}},
			{f: -0.1038, c: [4]int{2, 1, 0, 0}},
			{f: -0.0192, c: [4]int{2, 1, 0, -2}},
			{f: 0.0324, c: [4]int{2, 1, 0, -4}},
			{f: 0.0213, c: [4]int{2, -1, 0, 2}},
			{f: 0.1268, c: [4]int{2, -1, 0, 0}},
			{f: -0.0106, c: [4]int{1, 2, 0, 0}},
			{f: 0.0484, c: [4]int{1, 2, 0, -2}},
			{f: 0.0112, c: [4]int{1, -2, 0, 2}},
			{f: 0.0196, c: [4]int{1, -2, 0, 0}},
			{f: -0.0212, c: [4]int{1, -2, 0, -2}},
			{f: -0.0833, c: [4]int{1, 0, 2, -2}},
			{f: -0.0481, c: [4]int{1, 0, -2, 2}},
			{f: -0.7136, c: [4]int{1, 0, -2, 0}},
			{f: -0.0112, c: [4]int{1, 0, -2, -2}},
			{f: -0.01, c: [4]int{2, 0, 0, 1}},
			{f: 0.0155, c: [4]int{2, 0, 0, -1}},
			{f: 0.0164, c: [4]int{1, 1, 0, 1}},
			{f: 0.0401, c: [4]int{4, 0, 0, 0}},
			{f: -0.013, c: [4]int{4, 0, 0, -2}},
			{f: 0.0115, c: [4]int{3, -1, 0, 0}},
			{f: -0.0141, c: [4]int{2, 0, -2, -2}},
		},
	}
	sunf = [...][]float64{
		{
			-0.265, 0,
			3.76, 0,
			0.2, 0,
		},
		{
			-0.021, 0,
			5.18, 0,
			1.882, 3.8991,
			-0.03, 0,
		},
		{
			0.075, 5.1766,
			4.838, 5.2203,
			0.074, 3.6285,
			0.116, 2.5988,
			5.526, 2.5885,
			2.497, 5.5143,
			0.044, 5.435,
			0.666, 3.1016,
			1.559, 6.0258,
			1.024, 5.5527,
			0.21, 3.5989,
			0.144, 3.4104,
			0.152, 6.0004,
			0.084, 4.112,
			0.037, 3.8711,
			0.123, 3.4086,
			0.154, 6.2762,
			0.038, 4.6094,
			0.02, 5.1313,
			0.042, 4.5239,
			0.032, 0.8517,
			0.273, 3.7996,
			0.048, 4.5431,
			0.041, 6.0388,
			2.043, 6.002,
			1.77, 3.4977,
			0.028, 2.5831,
			0.129, 5.1348,
			0.425, 5.9146,
			0.034, 1.2391,
			0.5, 1.8357,
			0.585, 5.8304,
			0.085, 0.9529,
			0.204, 1.7593,
			0.02, 3.2463,
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
			2.6, 4.594,
			0.073, 4.8223,
			0.069, 1.4102,
			2.731, 1.521,
			1.61, 1.911,
			0.073, 4.4087,
			0.164, 2.9758,
			0.556, 1.4425,
			0.21, 1.7261,
			0.044, 2.9356,
			0.08, 1.3561,
			0.419, 1.7555,
			0.32, 4.703,
			0.108, 5.0719,
			0.112, 5.1243,
			0.021, 5.044,
		},
		{
			6.454, 0,
			0.177, 0,
			-0.424, 0,
			0.039, 0,
			-0.064, 0,
			0.172, 0,
		},
		{
			-0.092, 1.6354,
			-0.067, 2.1468,
			-0.21, 2.6494,
			-0.166, 4.6338,
		},
		{
			0.576, 0,
			-0.047, 0,
			0.021, 0,
		},
		{
			2.359e-6, 3.6607,
			6.842e-6, 1.018,
			0.869e-6, 3.9567,
			1.045e-6, 1.5332,
			1.497e-6, 4.4691,
			0.376e-6, 2.0295,
			2.057e-6, 4.0941,
			0.215e-6, 4.3459,
			0.478e-6, 0.2648,
			0.208e-6, 1.9548,
			7.067e-6, 1.563,
			0.244e-6, 5.9097,
			4.026e-6, 6.2526,
			1.459e-6, 0.3409,
			0.281e-6, 1.4172,
			0.803e-6, 6.1533,
			0.429e-6, 0.185,
		},
		{
			13.36e-6, 0,
			-1.33e-6, 0,
			0.37e-6, 0,
			0.36e-6, 0,
		},
	}
	sunc = [...][]int{
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
	mercf = [...][]float64{
		{
			0.013, 0.6807,
			0.048, 0.6283,
			0.185, 0.6231,
			0.711, 0.6191,
			0.285, 0.5784,
			0.075, 0.5411,
			0.019, 0.5585,
			0.01, 2.8449,
			0.039, 2.8117,
			0.147, 2.8135,
			0.552, 2.8126,
			2.1, 2.8126,
			3.724, 2.8046,
			0.729, 2.7883,
			0.186, 2.789,
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
			0.294, 3.991,
			0.484, 3.9514,
			0.07, 3.927,
			0.018, 3.927,
			0.013, 6.1261,
			0.05, 6.1052,
			0.185, 6.1069,
			0.685, 6.1011,
			2.81, 6.1062,
			7.356, 6.0699,
			1.471, 6.0685,
			0.375, 6.0687,
			0.098, 6.072,
			0.026, 6.0476,
			0.062, 5.154,
			0.122, 5.1191,
			0.011, 0.9076,
			0.074, 1.0123,
			0.106, 0.9372,
			0.017, 0.9425,
			0.02, 0.0506,
			0.052, 0.0384,
			0.052, 3.0281,
			0.012, 3.0543,
			0.011, 2.1642,
			0.016, 2.234,
			0.04, 4.3912,
			0.08, 4.4262,
			0.016, 4.4506,
		},
		{
			0.014, 1.0996,
			0.056, 1.1153,
			0.219, 1.116,
			0.083, 1.0734,
			0.024, 0.9442,
			0.018, 3.8432,
			0.07, 3.8293,
			0.256, 3.823,
			0.443, 3.8132,
			0.08, 3.7647,
			0.02, 3.7734,
			0.019, 0,
			0.133, 0.1134,
			0.129, 6.2588,
			0.026, 6.2413,
			0.026, 2.6599,
			0.087, 2.6232,
			0.374, 2.6496,
			0.808, 2.547,
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
			0.018, 0.665,
			0.069, 0.6405,
			0.253, 0.6449,
			0.938, 0.6454,
			3.275, 0.6458,
			0.499, 0.5569,
			0.119, 0.5271,
			0.032, 0.5184,
			0.03, 0.4939,
			0.106, 0.4171,
			0.353, 0.451,
			0.056, 0.384,
			0.013, 0.3142,
			0.028, 0.2531,
		},
		{
			0.034, 0.9512,
			0.06, 4.7962,
			0.028, 4.7124,
			0.028, 4.1836,
			0.102, 4.1871,
			0.38, 4.1864,
			0.059, 4.1818,
			0.015, 4.2185,
			0.012, 4.1713,
			0.05, 4.187,
		},
		{
			0.218e-6, 5.3369,
			0.491e-6, 5.3281,
			0.172e-6, 2.1642,
			0.091e-6, 2.1084,
			0.204e-6, 1.246,
			0.712e-6, 1.2413,
			2.370e-6, 1.2425,
			0.899e-6, 1.2303,
			0.763e-6, 4.3633,
			0.236e-6, 4.359,
			0.163e-6, 0.2705,
			0.541e-6, 0.271,
			1.157e-6, 0.259,
			0.099e-6, 0.1798,
			0.360e-6, 2.4237,
			0.234e-6, 2.374,
			0.253e-6, 4.5365,
			0.849e-6, 4.5293,
			2.954e-6, 4.5364,
			0.282e-6, 4.4581,
			1.550e-6, 1.357,
			0.472e-6, 1.3561,
			0.135e-6, 1.3579,
			0.081e-6, 3.5936,
			0.087e-6, 3.55,
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
			0.143e-6, 4.098,
		},
		{
			0.222e-6, 1.6024,
			0.708e-6, 1.5949,
			0.191e-6, 5.7914,
			0.1e-6, 5.3564,
			0.347e-6, 5.3548,
			1.185e-6, 5.3576,
			3.268e-6, 5.3579,
			0.371e-6, 2.2148,
			0.16e-6, 2.1241,
			0.134e-6, 5.126,
			0.347e-6, 5.162,
		},
	}
	mercc = [...][]int{
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
	venf = [...][]float64{
		{
			4.889, 2.0788,
			11.261, 2.587,
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
			0.3, 0.0218,
			0.159, 1.3491,
		},
		{
			2.246e-6, 0.508,
			9.772e-6, 1.0159,
			8.271e-6, 4.6674,
			0.737e-6, 0.8267,
			1.426e-6, 5.1747,
			0.51e-6, 5.7009,
			1.572e-6, 1.8188,
			0.717e-6, 2.2969,
			2.991e-6, 2.0611,
			1.335e-6, 0.9628,
		},
	}
	venc = [...][]int{
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
	nutf = [...][]float64{
		{
			0.2088, 0,
			-1.273, 0,
			0.1258, 0,
			-0.0496, 0,
			0.0214, 0,
			0.0124, 0,
		},
		{
			9.2109, 0,
			-0.0904, 0,
			0.5519, 0,
			0.0215, 0,
			-0.0093, 0,
			-0.0066, 0,
		},
		{
			-0.2037, 0,
			0.0675, 0,
			-0.0342, 0,
			-0.0261, 0,
			-0.0149, 0,
			0.0114, 0,
			0.006, 0,
			0.0058, 0,
			-0.0057, 0,
			-0.0052, 0,
		},
		{
			0.0884, 0,
			0.0183, 0,
			0.0113, 0,
			-0.005, 0,
		},
	}
	nutc = [...][]int{
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
		36525,       // eday of epoc
		19.19126393, // semi major axis (au)
		0.04716771,  // eccentricity
		0.76986,     // inclination (deg)
		74.22988,    // longitude of the ascending node (deg)
		170.96424,   // longitude of perihelion (deg)
		313.23218,   // mean longitude (deg)
		0.00152025,  // (au/julian century)
		-0.0001915,  // (e/julian century)
		-2.09,       // (arcsec/julian century)
		-1681.4,     // (arcsec/julian century)
		1312.56,     // (arcsec/julian century)
		1542547.79,  // (arcsec/julian century)
	}
	elemNept = []float64{
		36525,       // eday of epoc
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
		36525,       // eday of epoc
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
		522747.9,    // (arcsec/julian century)
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
	oStar = obj2{name: "Star", f: star}
	eObjs := [2]*obj2{{}, {}}
	if *eclipse != "" {
		if _, err := fmt.Sscanf(*eclipse, "%s %s", &eObjs[0].name, &eObjs[1].name); err != nil {
			log.Fatal("failed to parse eclipse objects")
		}
		for i, e := range eObjs {
			j := slices.IndexFunc(objs, func(o *obj2) bool {
				return e.name == o.name || e.name == o.fname
			})
			if j < 0 {
				log.Fatal("failed to parse eclipse objects")
			}
			eObjs[i] = objs[j]
		}
	}
	if root == "" {
		root = "/usr/local/plan9"
	}
	// Murray Hill, NJ.
	nlat = (40 + 41.06/60) * radian
	wlong = (74 + 23.98/60) * radian
	elev = 150 * metersToFeet
	if *loc != "" {
		if err = parseLocation(*loc); err != nil {
			log.Fatal("failed to parse location")
		}
	} else {
		data, err := os.ReadFile(filepath.Join(root, "sky", "here"))
		if err == nil {
			_ = parseLocation(string(data))
		}
	}
	glat = nlat - (692.74*radsec)*math.Sin(2*nlat) + (1.16*radsec)*math.Sin(4*nlat)
	erad = 0.99832707e0 + 0.00167644e0*math.Cos(2*nlat) - 0.352e-5*math.Cos(4*nlat) + 0.001e-5*math.Cos(6*nlat) + 0.1568e-6*elev
	for range *periods {
		d := day
		fmt.Print(julianToTime(d))
		if *printPos || *eclipse != "" {
			psTime(d)
		}
		fmt.Println()
		for i := range objs[0].point {
			seTime(d)
			for j := range objs {
				objs[j].f()
				obj(&objs[j].point[i])
				if *printPos {
					if *includeComet && objs[j].name != "comet" {
						continue
					}
					output(objs[j].fname, objs[j].point[i])
				}
			}
			if *eclipse != "" {
				d = dist(eObjs[0].point[i], eObjs[1].point[i])
				fmt.Printf("dist %s to %s = %.4f\n", eObjs[0].fname, eObjs[1].fname, d)
			}
			if *printPos || *eclipse != "" {
				break
			}
			d += stepSize
		}
		if !*printPos && *eclipse == "" {
			if err := search(); err != nil {
				log.Fatal(err)
			}
		}
		day += *interval
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
	const meanSolarDaySeconds = 0.001704
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

func fSun() {
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
	argp := (281.220833 + 0.0000470684*eday + 0.000453*capt2 + 0.000003*capt3) * radian
	anom := 358.475845 + 0.9856002670*eday - 0.00015*capt2 - 0.000003*capt3
	motion = 0.9856473354
	dmoon := 350.737681 + 12.1907491914*eday - 0.001436*capt2
	gmoon := 11.250889 + 13.229350449*eday - 0.003212*capt2
	mmoon := 296.104608 + 13.0649924465*eday + 9.192e-3*capt2
	mven := (212.448 + 1.602121635*eday) * radian
	merth := (358.476 + 0.985600267*eday) * radian
	mmars := (319.59 + 0.524024095*eday) * radian
	mjup := (225.269 + 0.083082362*eday) * radian
	msat := (175.593 + 0.033450794*eday) * radian
	dmoon = math.Mod(dmoon, 360) * radian
	gmoon = math.Mod(gmoon, 360) * radian
	mmoon = math.Mod(mmoon, 360) * radian
	anom += cosAdd(sunf[0], sunc[0], []float64{mmars, merth, mven, mjup}) / 3600
	anom += sinAdd(sunf[1], sunc[1], []float64{mmars, merth, mven, mjup, 0.07884 * capt}) / 3600
	anom = math.Mod(anom, 360) * radian
	// Computation of elliptic orbit.
	lambda = anom + argp
	pturbl := (6910.057-17.24*capt-0.052*capt2)*math.Sin(anom) +
		(72.338-0.361*capt)*math.Sin(2.*anom) +
		(1.054-0.001*capt)*math.Sin(3*anom) + 0.018*math.Sin(4*anom)
	lambda += pturbl * radsec
	beta = 0
	lograd := (30.57e-6 - 0.15e-6*capt) -
		(7274.12e-6-18.14e-6*capt-0.05e-6*capt2)*math.Cos(anom) -
		(91.38e-6-0.46e-6*capt)*math.Cos(2*anom) -
		(1.45e-6-0.01e-6*capt)*math.Cos(3*anom) -
		0.02e-6*math.Cos(4*anom)
	pturbl = cosAdd(sunf[2], sunc[2], []float64{mmars, merth, mven, mjup, msat})
	pturbl += sinAdd(sunf[3], sunc[3], []float64{dmoon, mmoon, merth}) + 0.9
	pturbl *= radsec
	pturbb := cosAdd(sunf[4], sunc[4], []float64{merth, mven, mjup})
	pturbb += sinAdd(sunf[5], sunc[5], []float64{gmoon, mmoon, dmoon})
	pturbb *= radsec
	pturbr := cosAdd(sunf[6], sunc[6], []float64{mmars, merth, mven, mjup, msat})
	pturbr += cosAdd(sunf[7], sunc[7], []float64{dmoon, mmoon, merth})
	lambda += pturbl
	if lambda > twoPi {
		lambda -= twoPi
	}
	beta += pturbb
	lograd = (lograd + pturbr) * 2.30258509
	rad = 1 + lograd*(1+lograd*(0.5+lograd/6))
	motion *= radian / (rad * rad)
	semi = 961.182
	if *searchOccult {
		semi = 959.63
	}
	mag = -26.5
}

// helio converts from ecliptic heliocentric coordinates referred to the mean
// equinox of date to equatorial geocentric coordinates referred to the true
// equator and equinox.
func helio() {
	// Compute geocentric distance of object and light-time correction.
	xmp := rad * math.Cos(beta) * math.Cos(lambda)
	ymp := rad * math.Cos(beta) * math.Sin(lambda)
	zmp := rad * math.Sin(beta)
	rp := math.Sqrt((xmp+xms)*(xmp+xms) + (ymp+yms)*(ymp+yms) + (zmp+zms)*(zmp+zms))
	lmb2 = lambda - 0.0057756e0*rp*motion
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
	beta2 := math.Atan2(zmp, math.Sqrt(xmp*xmp+ymp*ymp))
	lmb2 += phi
	// Change to equatorial coordinates.
	xmp = rp * math.Cos(lmb2) * math.Cos(beta2)
	ymp = rp * (math.Sin(lmb2)*math.Cos(beta2)*math.Cos(tobliq) - math.Sin(tobliq)*math.Sin(beta2))
	zmp = rp * (math.Sin(lmb2)*math.Cos(beta2)*math.Sin(tobliq) + math.Cos(tobliq)*math.Sin(beta2))
	alpha = math.Atan2(ymp, xmp)
	delta = math.Atan2(zmp, math.Sqrt(xmp*xmp+ymp*ymp))
	hp = 8.794e0 * radsec / rp
	semi /= rp
	if rad > 0 && rad < 2e5 {
		mag += 2.17 * math.Log(rad*rp)
	}
}

// geo converts geocentric equatorial coordinates to topocentric equatorial and
// topocentric horizon coordinates. All are (usually) referred to the true equator.
func geo() {
	// Convert to local hour angle and declination.
	lha = gst - alpha - wlong
	// Compute diurnal parallax (requires geocentric latitude).
	sa := math.Cos(delta) * math.Sin(lha)
	ca := math.Cos(delta)*math.Cos(lha) - erad*math.Cos(glat)*math.Sin(hp)
	sd := math.Sin(delta) - erad*math.Sin(glat)*math.Sin(hp)
	lha = math.Atan2(sa, ca)
	decl2 = math.Atan2(sd, math.Sqrt(sa*sa+ca*ca))
	f := math.Sqrt(sa*sa + ca*ca + sd*sd)
	semi2 = semi / f
	ra = gst - lha - wlong
	ra = piNorm(ra)
	// Convert to horizon coordinates.
	sel := math.Sin(nlat)*math.Sin(decl2) + math.Cos(nlat)*math.Cos(decl2)*math.Cos(lha)
	el = math.Atan2(sel, pyth(sel)) / radian
	saz := math.Sin(lha) * math.Cos(decl2)
	caz := math.Cos(nlat)*math.Sin(decl2) - math.Sin(nlat)*math.Cos(decl2)*math.Cos(lha)
	az = (math.Pi + math.Atan2(saz, -caz)) / radian
}

func piNorm(a float64) float64 {
	return math.Mod(a+math.Pi, 2*math.Pi) - math.Pi
}

func pyth(x float64) float64 {
	return math.Sqrt(1 - min(x*x, 1))
}

var k1, k2, k3, k4, mnom, msun, noded, dmoon float64

func moon() {
	// The fundamental elements - all referred to the epoch of Jan 0.5, 1900 and
	// to the mean equinox of date.
	dlong := 270.434164 + 13.1763965268*eday - 0.001133*capt2 + 2e-6*capt3
	argp := 334.329556 + 0.1114040803*eday - 0.010325*capt2 - 12e-6*capt3
	node := 259.183275 - 0.0529539222*eday + 0.002078*capt2 + 2e-6*capt3
	lsun := 279.696678 + 0.9856473354*eday + 0.000303*capt2
	psun := 281.220833 + 0.0000470684*eday + 0.000453*capt2 + 3e-6*capt3
	dlong = math.Mod(dlong, 360)
	argp = math.Mod(argp, 360)
	node = math.Remainder(node, 360)
	lsun = math.Mod(lsun, 360)
	psun = math.Mod(psun, 360)
	eccm := 22639.55
	eccs := 0.01675104 - 0.0000418*capt
	cpe := 124.986
	chp := 3422.451
	// Some subsidiary elements - they are all longitudes and they are referred
	// to the epoch 1/0.5 1900 and to the fixed mean equinox of 1850.0.
	v0 := 342.069128 + 1.602130482*eday
	t0 := 98.998753 + 0.9856091138*eday
	m0 := 293.049675 + 0.5240329445*eday
	j0 := 237.352319 + 0.0830912295*eday
	// The following are periodic corrections to the fundamental elements and constants.
	arg1 := (41.1 + 20.2*(capt+0.5)) * radian
	arg2 := (dlong - argp + 33 + 3*t0 - 10*v0 - 2.6*(capt+0.5)) * radian
	arg3 := (dlong - argp + 151.1 + 16*t0 - 18*v0 - (capt + 0.5)) * radian // The "Great Venus Inequality".
	arg4 := node * radian
	arg5 := (node + 276.2 - 2.3*(capt+0.5)) * radian
	arg6 := (313.9 + 13*t0 - 8*v0) * radian
	arg7 := (dlong - argp + 112 + 29*t0 - 26*v0) * radian
	arg8 := (dlong + argp - 2*lsun + 273 + 21*t0 - 20*v0) * radian
	arg9 := (node + 290.1 - 0.9*(capt+0.5)) * radian
	arg10 := (115 + 38.5*(capt+0.5)) * radian
	dlong += (0.84*math.Sin(arg1) + 0.31*math.Sin(arg2) + 14.27*math.Sin(arg3) + 7.261*math.Sin(arg4) +
		0.282*math.Sin(arg5) + 0.237*math.Sin(arg6) + 0.108*math.Sin(arg7) + 0.126*math.Sin(arg8)) / 3600
	argp += (-2.1*math.Sin(arg1) - 0.118*math.Sin(arg3) - 2.076*math.Sin(arg4) - 0.84*math.Sin(arg5) - 0.593*math.Sin(arg6)) / 3600
	node += (0.63*math.Sin(arg1) + 0.17*math.Sin(arg3) + 95.96*math.Sin(arg4) + 15.58*math.Sin(arg5) + 1.86*math.Sin(arg9)) / 3600
	t0 += (-6.4*math.Sin(arg1) - 1.89*math.Sin(arg6)) / 3600
	psun += (6.4*math.Sin(arg1) + 1.89*math.Sin(arg6)) / 3600
	dgamma := -4.318*math.Cos(arg4) - 0.698*math.Cos(arg5) - 0.083*math.Cos(arg9)
	j0 += 0.33 * math.Sin(arg10)
	// The following factors account for the fact that the eccentricity, solar eccentricity,
	// inclination and parallax used by Brown to make up his coefficients are
	// both wrong and out of date. Brown did the same thing in a different way.
	k1 = eccm / 22639.5
	k2 = eccs / 0.01675104
	k3 = 1 + 2.708e-6 + 0.000108008*dgamma
	k4 = cpe / 125.154
	k5 := chp / 3422.7
	// The principal arguments that are used to compute perturbations are the following
	// differences of the fundamental elements.
	mnom = dlong - argp
	msun = lsun - psun
	noded = dlong - node
	dmoon = dlong - lsun
	// Solar terms in longitude.
	var lterms float64
	for _, t := range moonTabs[0] {
		lterms += sinx(t.f, t.c[0], t.c[1], t.c[2], t.c[3], 0)
	}
	// Planetary terms in longitude.
	lterms += sinx(0.822, 0, 0, 0, 0, t0-v0)
	lterms += sinx(0.307, 0, 0, 0, 0, 2*t0-2*v0+179.8)
	lterms += sinx(0.348, 0, 0, 0, 0, 3*t0-2*v0+272.9)
	lterms += sinx(0.176, 0, 0, 0, 0, 4*t0-3*v0+271.7)
	lterms += sinx(0.092, 0, 0, 0, 0, 5*t0-3*v0+199)
	lterms += sinx(0.129, 1, 0, 0, 0, -t0+v0+180)
	lterms += sinx(0.152, 1, 0, 0, 0, t0-v0)
	lterms += sinx(0.127, 1, 0, 0, 0, 3*t0-3*v0+180)
	lterms += sinx(0.099, 0, 0, 0, 2, t0-v0)
	lterms += sinx(0.136, 0, 0, 0, 2, 2*t0-2*v0+179.5)
	lterms += sinx(0.083, -1, 0, 0, 2, -4*t0+4*v0+180)
	lterms += sinx(0.662, -1, 0, 0, 2, -3*t0+3*v0+180)
	lterms += sinx(0.137, -1, 0, 0, 2, -2*t0+2*v0)
	lterms += sinx(0.133, -1, 0, 0, 2, t0-v0)
	lterms += sinx(0.157, -1, 0, 0, 2, 2*t0-2*v0+179.6)
	lterms += sinx(0.079, -1, 0, 0, 2, -8*t0+6*v0+162.6)
	lterms += sinx(0.073, 2, 0, 0, -2, 3*t0-3*v0+180)
	lterms += sinx(0.643, 0, 0, 0, 0, -t0+j0+178.8)
	lterms += sinx(0.187, 0, 0, 0, 0, -2*t0+2*j0+359.6)
	lterms += sinx(0.087, 0, 0, 0, 0, j0+289.9)
	lterms += sinx(0.165, 0, 0, 0, 0, -t0+2*j0+241.5)
	lterms += sinx(0.144, 1, 0, 0, 0, t0-j0+1)
	lterms += sinx(0.158, 1, 0, 0, 0, -t0+j0+179)
	lterms += sinx(0.19, 1, 0, 0, 0, -2*t0+2*j0+180)
	lterms += sinx(0.096, 1, 0, 0, 0, -2*t0+3*j0+352.5)
	lterms += sinx(0.07, 0, 0, 0, 2, 2*t0-2*j0+180)
	lterms += sinx(0.167, 0, 0, 0, 2, -t0+j0+178.5)
	lterms += sinx(0.085, 0, 0, 0, 2, -2*t0+2*j0+359.2)
	lterms += sinx(1.137, -1, 0, 0, 2, 2*t0-2*j0+180.3)
	lterms += sinx(0.211, -1, 0, 0, 2, -t0+j0+178.4)
	lterms += sinx(0.089, -1, 0, 0, 2, -2*t0+2*j0+359.2)
	lterms += sinx(0.436, -1, 0, 0, 2, 2*t0-3*j0+7.5)
	lterms += sinx(0.24, 2, 0, 0, -2, -2*t0+2*j0+179.9)
	lterms += sinx(0.284, 2, 0, 0, -2, -2*t0+3*j0+172.5)
	lterms += sinx(0.195, 0, 0, 0, 0, -2*t0+2*m0+180.2)
	lterms += sinx(0.327, 0, 0, 0, 0, -t0+2*m0+224.4)
	lterms += sinx(0.093, 0, 0, 0, 0, -2*t0+4*m0+244.8)
	lterms += sinx(0.073, 1, 0, 0, 0, -t0+2*m0+223.3)
	lterms += sinx(0.074, 1, 0, 0, 0, t0-2*m0+306.3)
	lterms += sinx(0.189, 0, 0, 0, 0, node+180)
	// Solar terms in latitude.
	var sterms float64
	for _, t := range moonTabs[1] {
		sterms += sinx(t.f, t.c[0], t.c[1], t.c[2], t.c[3], 0)
	}
	var cterms float64
	for _, t := range moonTabs[2] {
		cterms += cosx(t.f, t.c[0], t.c[1], t.c[2], t.c[3], 0)
	}
	var nterms float64
	for _, t := range moonTabs[3] {
		nterms += sinx(t.f, t.c[0], t.c[1], t.c[2], t.c[3], 0)
	}
	// Planetary terms in latitude.
	pterms := sinx(0.215, 0, 0, 0, 0, dlong)
	// Solar terms in parallax.
	spterms := 3422.7
	for _, t := range moonTabs[4] {
		spterms += cosx(t.f, t.c[0], t.c[1], t.c[2], t.c[3], 0)
	}
	// Computation of longitude.
	lambda = (dlong + lterms/3600) * radian
	// Computation of latitude.
	arglat := (noded + sterms/3600) * radian
	gamma1 := 18519.7 * k3
	gamma2 := -6.241 * k3 * k3 * k3
	gamma3 := 0.004 * k3 * k3 * k3 * k3 * k3
	k6 := (gamma1 + cterms) / gamma1
	beta = k6*(gamma1*math.Sin(arglat)+gamma2*math.Sin(3.*arglat)+gamma3*math.Sin(5.*arglat)+nterms) + pterms
	if *searchOccult {
		beta -= 0.6
	}
	beta *= radsec
	// Computation of parallax.
	spterms = k5 * spterms * radsec
	hp = spterms + (spterms*spterms*spterms)/6
	rad = hp / radsec
	semi = 0.0799 + 0.272453*(hp/radsec)
	if dmoon < 0 {
		dmoon += 360
	}
	mag = dmoon / 360
	// Change to equatorial coordinates.
	lambda += phi
	obl2 := obliq + eps
	xmp := math.Cos(lambda) * math.Cos(beta)
	ymp := math.Sin(lambda)*math.Cos(beta)*math.Cos(obl2) - math.Sin(obl2)*math.Sin(beta)
	zmp := math.Sin(lambda)*math.Cos(beta)*math.Sin(obl2) + math.Cos(obl2)*math.Sin(beta)
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
	if m&1 > 0 {
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
		fSun()
	}
	if meday != eday {
		moon()
	}
	alpha = math.Mod(salph+math.Pi, twoPi)
	delta = -sdelt
	hp = mhp
	semi = 1.0183*mhp/radsec - 969.85/srad
	geo()
}

func merc() {
	ecc := 0.20561421 + 0.00002046*capt - 0.03e-6*capt2
	incl := (7.0028806 + 0.0018608*capt - 18.3e-6*capt2) * radian
	node := (47.145944 + 1.185208*capt + 0.0001739*capt2) * radian
	argp := (75.899697 + 1.55549*capt + 0.0002947*capt2) * radian
	mrad := 0.3870986
	anom := 102.279381 + 4.0923344364*eday + 6.7e-6*capt2
	motion = 4.0923770233
	q0 := (102.28 + 4.092334429*eday) * radian
	v0 := (212.536 + 1.602126105*eday) * radian
	t0 := (-1.45 + 0.985604737*eday) * radian
	j0 := (225.36 + 0.083086735*eday) * radian
	s0 := (175.68 + 0.033455441*eday) * radian
	anom = math.Mod(anom, 360) * radian
	enom := anom + ecc*math.Sin(anom)
	for {
		dele := (anom - enom + ecc*math.Sin(enom)) / (1 - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom := 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	pturbl := cosAdd(mercf[0], mercc[0], []float64{q0, -v0})
	pturbl += cosAdd(mercf[1], mercc[1], []float64{q0, -t0})
	pturbl += cosAdd(mercf[2], mercc[2], []float64{q0, -j0})
	pturbl += cosAdd(mercf[3], mercc[3], []float64{q0, -s0})
	pturbr := cosAdd(mercf[4], mercc[4], []float64{q0, -v0})
	pturbr += cosAdd(mercf[5], mercc[5], []float64{q0, -t0})
	pturbr += cosAdd(mercf[6], mercc[6], []float64{q0, -j0})
	// Reduce to the ecliptic.
	lambda = vnom + argp + pturbl*radsec
	nd := lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl := math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, pyth(sl))
	lograd := pturbr * 2.30258509
	rad *= 1 + lograd
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 3.34
	lsun := (99.696678 + 0.9856473354*eday) * radian
	elong := lambda - lsun
	ci := (rad - math.Cos(elong)) / math.Sqrt(1+rad*rad-2*rad*math.Cos(elong))
	dlong := math.Atan2(pyth(ci), ci) / radian
	mag = -0.003 + 0.01815*dlong + 0.0001023*dlong*dlong
	helio()
	geo()
}

func cosAdd(caf []float64, cac []int, coefs []float64) float64 {
	return trigAdd(math.Cos, caf, cac, coefs)
}

func sinAdd(caf []float64, cac []int, coefs []float64) float64 {
	return trigAdd(math.Sin, caf, cac, coefs)
}

func trigAdd(f func(float64) float64, caf []float64, cac []int, coefs []float64) float64 {
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
	// Mean orbital elements.
	ecc := 0.00682069 - 0.00004774*capt + 0.091e-6*capt2
	incl := (3.393631 + 0.0010058*capt - 0.97e-6*capt2) * radian
	node := (75.779647 + 0.89985*capt + 0.00041*capt2) * radian
	argp := (130.163833 + 1.408036*capt - 0.0009763*capt2) * radian
	mrad := 0.7233316
	anom := 212.603219 + 1.602130154*eday + 0.00128605*capt2
	motion = 1.6021687039
	// Mean anomalies of perturbing planets.
	v0 := (212.6 + 1.602130154*eday) * radian
	t0 := (358.63 + 0.985608747*eday) * radian
	m0 := (319.74 + 0.52403249*eday) * radian
	j0 := (225.43 + 0.083090842*eday) * radian
	s0 := (175.8 + 0.033459258*eday) * radian
	anom = math.Mod(anom, 360) * radian
	// Computation of long period terms affecting the mean anomaly.
	anom += (2.761-0.022*capt)*radsec*math.Sin(13*t0-8*v0+43.83*radian+4.52*radian*capt) +
		0.268*radsec*math.Cos(4*m0-7*t0+3*v0) +
		0.019*radsec*math.Sin(4*m0-7*t0+3*v0) -
		0.208*radsec*math.Sin(s0+1.4*radian*capt)
		// Computation of elliptic orbit.
	enom := anom + ecc*math.Sin(anom)
	for {
		dele := (anom - enom + ecc*math.Sin(enom)) / (1 - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom := 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	lambda = vnom + argp
	// Perturbations in longitude.
	pturbl := cosAdd(venf[0], venc[0], []float64{v0, t0, m0, j0}) * radsec
	// Perturbations in latitude.
	pturbb := cosAdd(venf[1], venc[1], []float64{v0, t0, j0}) * radsec
	// Perturbations in log radius vector.
	pturbr := cosAdd(venf[2], venc[2], []float64{v0, t0, m0, j0})
	// Reduction to the ecliptic.
	lambda += pturbl
	nd := lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl := math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, pyth(sl)) + pturbb
	lograd := pturbr * 2.30258509
	rad *= 1 + lograd
	motion *= radian * mrad * mrad / (rad * rad)
	// Computation of magnitude.
	lsun := (99.696678 + 0.9856473354*eday) * radian
	elong := lambda - lsun
	ci := (rad - math.Cos(elong)) / math.Sqrt(1+rad*rad-2*rad*math.Cos(elong))
	dlong := math.Atan2(pyth(ci), ci) / radian
	mag = -4 + 0.01322*dlong + 0.0000004247*dlong*dlong*dlong
	semi = 8.41
	helio()
	geo()
}

func mars() {
	ecc := 0.0933129 + 0.000092064*capt
	incl := (1.850333 - 6.75e-4*capt) * radian
	node := (48.786442 + 0.770992*capt) * radian
	argp := (334.218203 + 1.840758*capt + 1.3e-4*capt2) * radian
	mrad := 1.5236915
	anom := 319.529425 + 0.5240207666*eday + 1.808e-4*capt2
	motion = 0.5240711638
	anom = math.Mod(anom, 360) * radian
	enom := anom + ecc*math.Sin(anom)
	for {
		dele := (anom - enom + ecc*math.Sin(enom)) / (1 - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom := 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	lambda = vnom + argp
	// Reduce to the ecliptic.
	nd := lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl := math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, pyth(sl))
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 4.68
	lsun := (99.696678 + 0.9856473354*eday) * radian
	elong := lambda - lsun
	ci := (rad - math.Cos(elong)) / math.Sqrt(1+rad*rad-2*rad*math.Cos(elong))
	dlong := math.Atan2(pyth(ci), ci) / radian
	mag = -1.3 + 0.01486*dlong
	helio()
	geo()
}

func jup() {
	ecc := 0.0483376 + 163e-6*capt
	incl := (1.30866 - 0.0055*capt) * radian
	node := (99.43785 + 1.011*capt) * radian
	argp := (12.71165 + 1.611*capt) * radian
	mrad := 5.202803
	anom := 225.22165 + 0.0830912*eday - 0.0484*capt
	motion = 299.1284 / 3600
	anom = math.Mod(anom, 360) * radian
	enom := anom + ecc*math.Sin(anom)
	for {
		dele := (anom - enom + ecc*math.Sin(enom)) / (1 - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom := 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	lambda = vnom + argp
	// Reduce to the ecliptic.
	nd := lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd)) + 555*radsec
	sl := math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, pyth(sl)) - 51*radsec
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 98.47
	mag = -8.93
	helio()
	geo()
}

func sat() {
	ecc := 0.05589 - 0.000347*capt
	incl := (2.49256 - 0.0044*capt) * radian
	node := (112.78364 + 0.87306*capt) * radian
	argp := (91.08897 + 1.95917*capt) * radian
	mrad := 9.538843
	anom := 175.4763 + 0.03345972*eday - 0.56527*capt
	motion = 120.455 / 3600
	anom = math.Mod(anom, 360) * radian
	enom := anom + ecc*math.Sin(anom)
	for {
		dele := (anom - enom + ecc*math.Sin(enom)) / (1 - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom := 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	lambda = vnom + argp
	// Reduce to the ecliptic.
	nd := lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd)) - 1185*radsec
	sl := math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, pyth(sl)) - 51*radsec
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 83.33
	// Computation of magnitude; first, find the geocentric equatorial coordinates of Saturn.
	sd := rad*(math.Cos(beta)*math.Sin(lambda)*math.Sin(obliq)+math.Sin(beta)*math.Cos(obliq)) + zms
	sa := rad*(math.Cos(beta)*math.Sin(lambda)*math.Cos(obliq)-math.Sin(beta)*math.Sin(obliq)) + yms
	ca := rad*(math.Cos(beta)*math.Cos(lambda)) + xms
	alpha = math.Atan2(sa, ca)
	delta = math.Atan2(sd, math.Sqrt(sa*sa+ca*ca))
	// Here are the necessary elements of Saturn's rings cf. Exp. Supp. p. 363ff.
	capj := (6.9056 - 0.4322*capt) * radian
	capn := (126.3615 + 3.9894*capt + 0.2403*capt2) * radian
	eye := (28.0743 - 0.0128*capt) * radian
	comg := (168.1179 + 1.3936*capt) * radian
	omg := (42.9236 - 2.739*capt - 0.2344*capt2) * radian
	// Now find saturnicentric ring-plane coords of the earth.
	sb := math.Sin(capj)*math.Cos(delta)*math.Sin(alpha-capn) - math.Cos(capj)*math.Sin(delta)
	su := math.Cos(capj)*math.Cos(delta)*math.Sin(alpha-capn) + math.Sin(capj)*math.Sin(delta)
	cu := math.Cos(delta) * math.Cos(alpha-capn)
	u := math.Atan2(su, cu)
	b := math.Atan2(sb, math.Sqrt(su*su+cu*cu))
	// And then the saturnicentric ring-plane coords of the sun.
	su = math.Sin(eye)*math.Sin(beta) + math.Cos(eye)*math.Cos(beta)*math.Sin(lambda-comg)
	cu = math.Cos(beta) * math.Cos(lambda-comg)
	up := math.Atan2(su, cu)
	// At last, the magnitude.
	sb = math.Sin(b)
	mag = -8.68 + 2.52*math.Abs(up+omg-u) - 2.6*math.Abs(sb) + 1.25*(sb*sb)
	helio()
	geo()
}

func uran() {
	cy := (eday - elemUran[0]) / 36525 // Per julian century.
	mrad := elemUran[1] + elemUran[1+6]*cy
	ecc := elemUran[2] + elemUran[2+6]*cy
	cy = cy / 3600 // arcsec/deg per julian century.
	incl := (elemUran[3] + elemUran[3+6]*cy) * radian
	node := (elemUran[4] + elemUran[4+6]*cy) * radian
	argp := (elemUran[5] + elemUran[5+6]*cy)
	anom := elemUran[6] + elemUran[6+6]*cy - argp
	motion = elemUran[6+6] / 36525 / 3600
	argp *= radian
	anom = math.Mod(anom, 360) * radian
	enom := anom + ecc*math.Sin(anom)
	for {
		dele := (anom - enom + ecc*math.Sin(enom)) / (1 - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom := 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	lambda = vnom + argp
	// Reduce to the ecliptic.
	nd := lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd)) - 1185*radsec
	sl := math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, pyth(sl)) - 51*radsec
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 83.33
	// Computation of magnitude; first, find the geocentric equatorial coordinates of Saturn.
	sd := rad*(math.Cos(beta)*math.Sin(lambda)*math.Sin(obliq)+math.Sin(beta)*math.Cos(obliq)) + zms
	sa := rad*(math.Cos(beta)*math.Sin(lambda)*math.Cos(obliq)-math.Sin(beta)*math.Sin(obliq)) + yms
	ca := rad*(math.Cos(beta)*math.Cos(lambda)) + xms
	alpha = math.Atan2(sa, ca)
	delta = math.Atan2(sd, math.Sqrt(sa*sa+ca*ca))
	// Here are the necessary elements of Saturn's rings cf. Exp. Supp. p. 363ff.
	capj := (6.9056 - 0.4322*capt) * radian
	capn := (126.3615 + 3.9894*capt + 0.2403*capt2) * radian
	eye := (28.0743 - 0.0128*capt) * radian
	comg := (168.1179 + 1.3936*capt) * radian
	omg := (42.9236 - 2.739*capt - 0.2344*capt2) * radian
	// Now find saturnicentric ring-plane coords of the earth.
	sb := math.Sin(capj)*math.Cos(delta)*math.Sin(alpha-capn) - math.Cos(capj)*math.Sin(delta)
	su := math.Cos(capj)*math.Cos(delta)*math.Sin(alpha-capn) + math.Sin(capj)*math.Sin(delta)
	cu := math.Cos(delta) * math.Cos(alpha-capn)
	u := math.Atan2(su, cu)
	b := math.Atan2(sb, math.Sqrt(su*su+cu*cu))
	// And then the saturnicentric ring-plane coords of the sun.
	su = math.Sin(eye)*math.Sin(beta) + math.Cos(eye)*math.Cos(beta)*math.Sin(lambda-comg)
	cu = math.Cos(beta) * math.Cos(lambda-comg)
	up := math.Atan2(su, cu)
	// At last, the magnitude.
	sb = math.Sin(b)
	mag = -8.68 + 2.52*math.Abs(up+omg-u) - 2.6*math.Abs(sb) + 1.25*(sb*sb)
	helio()
	geo()
}

func nept() {
	cy := (eday - elemNept[0]) / 36525 // Per julian century.
	mrad := elemNept[1] + elemNept[1+6]*cy
	ecc := elemNept[2] + elemNept[2+6]*cy
	cy = cy / 3600 // arcsec/deg per julian century.
	incl := (elemNept[3] + elemNept[3+6]*cy) * radian
	node := (elemNept[4] + elemNept[4+6]*cy) * radian
	argp := (elemNept[5] + elemNept[5+6]*cy)
	anom := elemNept[6] + elemNept[6+6]*cy - argp
	motion = elemNept[6+6] / 36525 / 3600
	argp *= radian
	anom = math.Mod(anom, 360) * radian
	enom := anom + ecc*math.Sin(anom)
	for {
		dele := (anom - enom + ecc*math.Sin(enom)) / (1 - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom := 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	lambda = vnom + argp
	// Reduce to the ecliptic.
	nd := lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl := math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, pyth(sl))
	lambda -= 1185 * radsec
	beta -= 51 * radsec
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 83.33
	// Computation of magnitude; first, find the geocentric equatorial coordinates of Saturn.
	sd := rad*(math.Cos(beta)*math.Sin(lambda)*math.Sin(obliq)+math.Sin(beta)*math.Cos(obliq)) + zms
	sa := rad*(math.Cos(beta)*math.Sin(lambda)*math.Cos(obliq)-math.Sin(beta)*math.Sin(obliq)) + yms
	ca := rad*(math.Cos(beta)*math.Cos(lambda)) + xms
	alpha = math.Atan2(sa, ca)
	delta = math.Atan2(sd, math.Sqrt(sa*sa+ca*ca))
	// Here are the necessary elements of Saturn's rings cf. Exp. Supp. p. 363ff.
	capj := (6.9056 - 0.4322*capt) * radian
	capn := (126.3615 + 3.9894*capt + 0.2403*capt2) * radian
	eye := (28.0743 - 0.0128*capt) * radian
	comg := (168.1179 + 1.3936*capt) * radian
	omg := (42.9236 - 2.739*capt - 0.2344*capt2) * radian
	// Now find saturnicentric ring-plane coords of the earth.
	sb := math.Sin(capj)*math.Cos(delta)*math.Sin(alpha-capn) - math.Cos(capj)*math.Sin(delta)
	su := math.Cos(capj)*math.Cos(delta)*math.Sin(alpha-capn) + math.Sin(capj)*math.Sin(delta)
	cu := math.Cos(delta) * math.Cos(alpha-capn)
	u := math.Atan2(su, cu)
	b := math.Atan2(sb, math.Sqrt(su*su+cu*cu))
	// And then the saturnicentric ring-plane coords of the sun.
	su = math.Sin(eye)*math.Sin(beta) + math.Cos(eye)*math.Cos(beta)*math.Sin(lambda-comg)
	cu = math.Cos(beta) * math.Cos(lambda-comg)
	up := math.Atan2(su, cu)
	// At last, the magnitude.
	sb = math.Sin(b)
	mag = -8.68 + 2.52*math.Abs(up+omg-u) - 2.6*math.Abs(sb) + 1.25*(sb*sb)
	helio()
	geo()
}

func plut() {
	cy := (eday - elemPlut[0]) / 36525 // Per julian century.
	mrad := elemPlut[1] + elemPlut[1+6]*cy
	ecc := elemPlut[2] + elemPlut[2+6]*cy
	cy = cy / 3600 // arcsec/deg per julian century.
	incl := (elemPlut[3] + elemPlut[3+6]*cy) * radian
	node := (elemPlut[4] + elemPlut[4+6]*cy) * radian
	argp := (elemPlut[5] + elemPlut[5+6]*cy)
	anom := elemPlut[6] + elemPlut[6+6]*cy - argp
	motion = elemPlut[6+6] / 36525 / 3600
	argp *= radian
	anom = math.Mod(anom, 360) * radian
	enom := anom + ecc*math.Sin(anom)
	for {
		dele := (anom - enom + ecc*math.Sin(enom)) / (1 - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom := 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	lambda = vnom + argp
	// Reduce to the ecliptic.
	nd := lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl := math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, pyth(sl))
	lambda -= 1185 * radsec
	beta -= 51 * radsec
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 83.33
	// Computation of magnitude; first, find the geocentric equatorial coordinates of Saturn.
	sd := rad*(math.Cos(beta)*math.Sin(lambda)*math.Sin(obliq)+math.Sin(beta)*math.Cos(obliq)) + zms
	sa := rad*(math.Cos(beta)*math.Sin(lambda)*math.Cos(obliq)-math.Sin(beta)*math.Sin(obliq)) + yms
	ca := rad*(math.Cos(beta)*math.Cos(lambda)) + xms
	alpha = math.Atan2(sa, ca)
	delta = math.Atan2(sd, math.Sqrt(sa*sa+ca*ca))
	// Here are the necessary elements of Saturn's rings cf. Exp. Supp. p. 363ff.
	capj := (6.9056 - 0.4322*capt) * radian
	capn := (126.3615 + 3.9894*capt + 0.2403*capt2) * radian
	eye := (28.0743 - 0.0128*capt) * radian
	comg := (168.1179 + 1.3936*capt) * radian
	omg := (42.9236 - 2.739*capt - 0.2344*capt2) * radian
	// Now find saturnicentric ring-plane coords of the earth.
	sb := math.Sin(capj)*math.Cos(delta)*math.Sin(alpha-capn) - math.Cos(capj)*math.Sin(delta)
	su := math.Cos(capj)*math.Cos(delta)*math.Sin(alpha-capn) + math.Sin(capj)*math.Sin(delta)
	cu := math.Cos(delta) * math.Cos(alpha-capn)
	u := math.Atan2(su, cu)
	b := math.Atan2(sb, math.Sqrt(su*su+cu*cu))
	// And then the saturnicentric ring-plane coords of the sun.
	su = math.Sin(eye)*math.Sin(beta) + math.Cos(eye)*math.Cos(beta)*math.Sin(lambda-comg)
	cu = math.Cos(beta) * math.Cos(lambda-comg)
	up := math.Atan2(su, cu)
	// At last, the magnitude.
	sb = math.Sin(b)
	mag = -8.68 + 2.52*math.Abs(up+omg-u) - 2.6*math.Abs(sb) + 1.25*(sb*sb)
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
	// 153P/Ikeya–Zhang.
	t := time.Date(2002, 3, 18, 23, 28, 53, 760000000, time.UTC)
	elem := cometElem{t: timeToJulian(&t) + 2415020, q: 0.5070601, e: 0.990111, i: 28.12106, w: 34.6666, o: 93.1206}
	ecc := elem.e
	if ecc > 0.999 { // Can't do hyperbolas.
		ecc = 0.999
	}
	incl := elem.i * radian
	node := (elem.o + 0.4593) * radian
	argp := (elem.w + elem.o + 0.4066) * radian
	mrad := elem.q / (1 - ecc)
	motion = 0.01720209895 * math.Sqrt(1/(mrad*mrad*mrad)) / radian
	anom := (eday - (elem.t - 2415020)) * motion * radian
	enom := anom + ecc*math.Sin(anom)
	for {
		dele := (anom - enom + ecc*math.Sin(enom)) / (1 - ecc*math.Cos(enom))
		enom += dele
		if math.Abs(dele) <= converge {
			break
		}
	}
	vnom := 2 * math.Atan2(math.Sqrt((1+ecc)/(1-ecc))*math.Sin(enom/2), math.Cos(enom/2))
	rad = mrad * (1 - ecc*math.Cos(enom))
	lambda = vnom + argp
	// Reduce to the ecliptic.
	nd := lambda - node
	lambda = node + math.Atan2(math.Sin(nd)*math.Cos(incl), math.Cos(nd))
	sl := math.Sin(incl) * math.Sin(nd)
	beta = math.Atan2(sl, math.Sqrt(1-sl*sl))
	motion *= radian * mrad * mrad / (rad * rad)
	semi = 0
	mag = 5.47 + 6.1/2.303*math.Log(rad)
	helio()
	geo()
}

func star() {
	ra = oStar.point[0].ra
	decl2 = oStar.point[0].decl2
	semi2 = oStar.point[0].semi2
	az = oStar.point[0].az
	el = oStar.point[0].el
	mag = oStar.point[0].mag
}

func psTime(d float64) {
	seTime(d)
	semi = 0
	motion = 0
	rad = 1e9
	lambda = 0
	beta = 0
	helio()
	geo()
	fmt.Printf(" %s %s %s %4.0f", rConv(lha), dConv(nlat), dConv(awlong), elev/metersToFeet)
}

func seTime(d float64) {
	eday = d + ΔT/86400
	wlong = awlong + 15*ΔT*radsec
	capt = eday / 36524.22e0
	capt2 = capt * capt
	capt3 = capt * capt2
	nutate()
	eday += 0.1
	sun()
	srad = rad
	xm := rad * math.Cos(beta) * math.Cos(lambda)
	ym := rad * math.Cos(beta) * math.Sin(lambda)
	zm := rad * math.Sin(beta)
	eday -= 0.1
	sun()
	xms = rad * math.Cos(beta) * math.Cos(lambda)
	yms = rad * math.Cos(beta) * math.Sin(lambda)
	zms = rad * math.Sin(beta)
	x := 0.057756
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
	mnom := (296.104608 + 13.0649924465*eday + 9.192e-3*capt2 + 14.38e-6*capt3) * radian
	msun := (358.475833 + 0.9856002669*eday - 0.15e-3*capt2 - 3.33e-6*capt3) * radian
	noded = (11.250889 + 13.229350449*eday - 3.211e-3*capt2 - 0.33e-6*capt3) * radian
	dmoon := (350.737486 + 12.1907491914*eday - 1.436e-3*capt2 + 1.89e-6*capt3) * radian
	node := (259.183275 - 0.0529539222*eday + 2.078e-3*capt2 + 2.22e-6*capt3) * radian
	phi = -(17.2327 + 0.01737*capt) * math.Sin(node)
	phi += sinAdd(nutf[0], nutc[0], []float64{node, noded, dmoon, msun})
	eps = cosAdd(nutf[1], nutc[1], []float64{node, noded, dmoon, msun})
	dphi := sinAdd(nutf[2], nutc[2], []float64{node, noded, mnom, dmoon})
	deps := cosAdd(nutf[3], nutc[3], []float64{node, noded, mnom})
	phi = (phi + dphi) * radsec
	eps = (eps + deps) * radsec
	obliq = (23.452294 - 0.0130125*capt - 1.64e-6*capt2 + 0.503e-6*capt3) * radian
	tobliq = obliq + eps
	gst = 99.690983 + 360.9856473354*eday + 0.000387*capt2 - 180
	gst = math.Mod(gst, 360)
	if gst < 0 {
		gst += 360
	}
	gst *= radian
	gst += phi * math.Cos(obliq)
}

func obj(o *obj1) {
	*o = obj1{ra: ra, decl2: decl2, semi2: semi2, az: az, el: el, mag: mag}
}

func output(n string, p obj1) {
	if n == "" {
		fmt.Printf(" SAO %s", sao)
	} else {
		fmt.Printf("%10s", n)
	}
	fmt.Printf(" %s %s %9.4f %9.4f %9.4f", rConv(p.ra), dConv(p.decl2), p.az, p.el, p.semi2)
	if n == oSun.fname || n == oMoon.fname {
		fmt.Printf(" %7.4f", p.mag)
	}
	fmt.Println()
}

func rConv(v float64) string {
	v = math.Mod(v*12/math.Pi+24, 24)
	h := math.Floor(v)
	v = math.Mod((v-h)*60, 60)
	m := math.Floor(v)
	v = math.Mod((v-m)*60, 60)
	c := math.Floor(v)
	return fmt.Sprintf("%2dh%.2dm%.2ds", int(h), int(m), int(c))
}

func dConv(v float64) string {
	v = math.Mod(v/radian, 360)
	if v < 0 {
		v += 360
	}
	sign := '+'
	if v > 180 {
		v = 360 - v
		sign = '-'
	}
	h := math.Floor(v)
	v = math.Mod((v-h)*60, 60)
	m := math.Floor(v)
	v = math.Mod((v-m)*60, 60)
	c := math.Floor(v)
	return fmt.Sprintf(`%c%.2d°%.2d'%.2d"`, sign, int(h), int(m), int(c))
}

func dist(o1, o2 obj1) float64 {
	d := math.Sin(o1.decl2)*math.Sin(o2.decl2) +
		math.Cos(o1.decl2)*math.Cos(o2.decl2)*math.Cos(o1.ra-o2.ra)
	return math.Abs(math.Atan2(pyth(d), d)) / radsec
}

func search() error {
	for i, o := range objs {
		if o.name == oShad.name {
			continue
		}
		t := rise(*o, -0.833)
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
		t = set(*o, -0.833)
		if t >= 0 {
			err := event(evt{s: fmt.Sprintf("%s sets at ", o.fname), tim: t, flag: flag})
			if err != nil {
				return err
			}
		}
		if o.name == oSun.name {
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
				{0.762, "Eta aquarid"},
				{1.5497, "Ophiuchid"},
				{2.1324, "Capricornid"},
				{2.1991, "Delta aquarid"},
				{2.2158, "Pisces australid"},
				{2.4331, "Perseid"},
				{-2.6578, "Orionid"},
				{-1.8678, "Phoenicid"},
				{-1.726, "Geminid"},
			}
			for _, b := range bettab {
				t = betCross(b.beta)
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
		if o.name == oMoon.name {
			for j := range len(o.point) - 2 {
				if o.point[j].mag > 0.75 && o.point[j+1].mag < 0.25 {
					err := event(evt{s: "New moon"})
					if err != nil {
						return err
					}
				}
				if o.point[j].mag <= 0.25 && o.point[j+1].mag > 0.25 {
					err := event(evt{s: "First quarter moon"})
					if err != nil {
						return err
					}
				}
				if o.point[j].mag <= 0.5 && o.point[j+1].mag > 0.5 {
					err := event(evt{s: "Full moon"})
					if err != nil {
						return err
					}
				}
				if o.point[j].mag <= 0.75 && o.point[j+1].mag > 0.75 {
					err := event(evt{s: "Last quarter moon"})
					if err != nil {
						return err
					}
				}
			}
		}
		if o.name == oMerc.name || o.name == oVenus.name {
			t = float64(meLong(*o))
			if t >= 0 {
				t = rise(*o, 0) - rise(oSun, 0)
				if t < 0 {
					t += numPoints
				}
				if t > numPoints {
					t -= numPoints
				}
				if t > numPoints/2 {
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
			if o.name == oMoon.name || p.name == oMoon.name {
				if err := occult(*o, *p); err != nil {
					return err
				}
				if occ.t3 < 0 {
					continue
				}
				if o.name == oSun.name || p.name == oMoon.name {
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
			if o.name == oSun.name {
				if p.name != oMerc.name && p.name != oVenus.name {
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
	if *searchOccult {
		if err := stars(); err != nil {
			return err
		}
	}
	flushEvents()
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
	for i := range len(oSun.point) - 1 {
		d1, d2 = d2, oSun.point[i].ra
		if n == 0 {
			d2 -= math.Pi
			if d2 < -math.Pi {
				d2 += twoPi
			}
		}
		if i >= 1 && d3 >= d1 && d3 < d2 {
			return float64(i) - (d3-d2)/(d1-d2)
		}
	}
	return -1
}

func betCross(b float64) float64 {
	for i := 1; i < len(oSun.point); i++ {
		d1, d2 := oSun.point[i-1].mag, oSun.point[i].mag
		if b >= d1 && b < d2 {
			return float64(i) - (b-d2)/(d1-d2)
		}
	}
	return -1
}

func meLong(o obj2) int {
	for i := 2; i < len(o.point); i++ {
		d1 := dist(o.point[i-2], oSun.point[i-2])
		d2 := dist(o.point[i-1], oSun.point[i-1])
		d3 := dist(o.point[i], oSun.point[i])
		if d2 >= d1 && d2 >= d3 {
			return i - 2
		}
	}
	return -1
}

func event(e evt) error {
	if e.flag&dark > 0 && sunEl(e.tim) > -12 {
		return nil
	}
	if e.flag&light > 0 && sunEl(e.tim) < 0 {
		return nil
	}
	if len(events) >= 100 {
		return errors.New("too many events")
	}
	events = append(events, e)
	return nil
}

func flushEvents() {
	slices.SortFunc(events, func(e1, e2 evt) int {
		t1, t2 := e1.tim, e2.tim
		if e1.flag&signif > 0 {
			t1 -= 1000
		}
		if e2.flag&signif > 0 {
			t2 -= 1000
		}
		return cmp.Compare(t1, t2)
	})
	for _, e := range events {
		if e.flag&ptime > 0 {
			fmt.Printf("%s%s\n", e.s, julianToTime(day+e.tim*stepSize).Format(time.TimeOnly+" MST"))
		} else {
			fmt.Printf("%s\n", e.s)
		}
	}
	events = nil
}

func sunEl(t float64) float64 {
	i := int(t)
	if i < 0 || i > numPoints {
		return -90
	}
	return oSun.point[i].el + (t-float64(i))*(oSun.point[i+1].el-oSun.point[i].el)
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
	n := 2880 * iVal / numPoints // 1 min steps.
	i -= 2
	pts(o1, i, &occ1)
	pts(o2, i, &occ2)
	di := float64(i)
	x := 0.
	dx := 2 / n
	for i = range int(n + 1) {
		pt(&occ1, x)
		pt(&occ2, x)
		d1, d2 = d2, d3
		d3 = dist(occ1.act, occ2.act)
		if i >= 2 && d2 <= d1 && d2 <= d3 {
			ok = true
			break
		}
		x += dx
	}
	if !ok {
		return errors.New("bad 1\n")
	}
	ok = false
	if d2 > occ1.act.semi2+occ2.act.semi2+50 {
		return nil
	}
	di += x - 3*dx
	x = 0
	var xo1, xo2 obj2
	for i = range 3 {
		seTime(day + stepSize*(di+x))
		o1.f()
		obj(&xo1.point[i])
		o2.f()
		obj(&xo2.point[i])
		x += 2 * dx
	}
	dx /= 60
	x = 0
	pts(xo1, 0, &occ1)
	pts(xo2, 0, &occ2)
	for i = range 241 {
		pt(&occ1, x)
		pt(&occ2, x)
		d1, d2 = d2, d3
		d3 = dist(occ1.act, occ2.act)
		if i >= 2 && d2 <= d1 && d2 <= d3 {
			ok = true
			break
		}
		x += 1. / 120
	}
	if !ok {
		return errors.New("bad 2\n")
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
		pt(&occ1, x)
		pt(&occ2, x)
		d1 = d2
		d2 = dist(occ1.act, occ2.act)
		if i != i1 {
			if d1 <= d3 && d2 > d3 {
				occ.t4 = di + (float64(i)-0.5)*dx
				occ.e4 = occ1.act.el
			}
			if d2 > d4 {
				if d1 <= d4 {
					occ.t5 = di + (float64(i)-0.5)*dx
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
		pt(&occ1, x)
		pt(&occ2, x)
		d1 = d2
		d2 = dist(occ1.act, occ2.act)
		if i != i1 {
			if d1 <= d3 && d2 > d3 {
				occ.t2 = di + (float64(i)-0.5)*dx
				occ.e2 = occ1.act.el
			}
			if d2 > d4 {
				if d1 <= d4 {
					occ.t1 = di + (float64(i)-0.5)*dx
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

func pts(o obj2, i int, oc *occt) {
	p1, p2, p3 := o.point[i], o.point[i+1], o.point[i+2]
	oc.del0.ra = p1.ra
	oc.del0.decl2 = p1.decl2
	oc.del0.semi2 = p1.semi2
	oc.del0.el = p1.el
	a := p2.ra - p1.ra
	oc.del1.ra = piNorm(a)
	a = p2.decl2 - p1.decl2
	oc.del1.decl2 = piNorm(a)
	oc.del1.semi2 = p2.semi2 - p1.semi2
	oc.del1.el = p2.el - p1.el
	a = p1.ra + p3.ra - 2*p2.ra
	oc.del2.ra = piNorm(a) / 2
	a = p1.decl2 + p3.decl2 - 2*p2.decl2
	oc.del2.decl2 = piNorm(a) / 2
	oc.del2.semi2 = (p1.semi2 + p3.semi2 - 2*p2.semi2) / 2
	oc.del2.el = (p1.el + p3.el - 2*p2.el) / 2
}

func pt(o *occt, x float64) {
	y := x * (x - 1)
	o.act.ra = o.del0.ra + x*o.del1.ra + y*o.del2.ra
	o.act.decl2 = o.del0.decl2 + x*o.del1.decl2 + y*o.del2.decl2
	o.act.semi2 = o.del0.semi2 + x*o.del1.semi2 + y*o.del2.semi2
	o.act.el = o.del0.el + x*o.del1.el + y*o.del2.el
}

func stars() error {
	sd := 1000 * radsec
	lomoon := oMoon.point[0].ra - sd
	if lomoon < 0 {
		lomoon += twoPi
	}
	himoon := oMoon.point[numPoints+1].ra + sd
	if himoon > twoPi {
		himoon -= twoPi
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
	epoch := (1950-1900)*365.2422 + 0.313
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
		if wrap && (alpha < lomoon && alpha > himoon) || !wrap && (alpha < lomoon || alpha > himoon) {
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
		mag, err := strconv.ParseFloat(strings.TrimSpace(l[63:67]), 64)
		if err != nil {
			return err
		}
		// Convert rt ascension and declination to internal format.
		delta = float64(abs(dday)) + float64(dmin)/60 + dsec/3600
		if dday < 0 {
			delta = -delta
		}
		// Remove E-terms of aberration except when finding catalog mean places.
		alpha += (0.341 / (3600 * 15)) * math.Sin((alpha+11.26)*15*radian) / math.Cos(delta*radian)
		delta += (0.341/3600)*math.Cos((alpha+11.26)*15*radian)*math.Sin(delta*radian) - (0.029/3600)*math.Cos(delta*radian)
		// Correct for proper motion.
		tau := (eday - epoch) / 365.2422
		alpha += tau * da / 3600
		delta += tau * dd / 3600
		alpha *= 15 * radian
		delta *= radian
		// Convert to rectangular coordinates merely for convenience.
		xm := math.Cos(delta) * math.Cos(alpha)
		ym := math.Cos(delta) * math.Sin(alpha)
		zm := math.Sin(delta)
		// Convert mean places at epoch of startable to current epoch (i.e. compute relevant precession).
		capt0 := (epoch - 18262.427) / 36524.22e0
		capt1 := (eday - epoch) / 36524.22
		capt12 := capt1 * capt1
		capt13 := capt12 * capt1
		xx := -(0.00029696+26.e-8*capt0)*capt12 - 13e-8*capt13
		yx := -(0.02234941+1355.e-8*capt0)*capt1 - 676e-8*capt12 + 221e-8*capt13
		zx := -(0.0097169-414e-8*capt0)*capt1 + 207e-8*capt12 + 96e-8*capt13
		yy := -(0.00024975+30e-8*capt0)*capt12 - 15e-8*capt13
		zy := -(0.00010858 + 2e-8*capt0) * capt12
		zz := -(0.00004721 - 4e-8*capt0) * capt12
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
		rad = 1e9
		if px != 0 {
			rad = 20600 / px
		}
		motion = 0
		semi = 0
		helio()
		geo()
		sd = 0.0896833e0*math.Cos(beta)*math.Sin(lambda-1.382+.00092422117*eday) + 0.99597*math.Sin(beta)
		if math.Abs(sd) > 0.0183 {
			continue
		}
		for i := range oStar.point {
			obj(&oStar.point[i])
		}
		if err = occult(oMoon, oStar); err != nil {
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
