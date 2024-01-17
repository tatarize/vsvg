//! Recreation of Georg Nees' ["Schotter" (1968-1970)](https://collections.vam.ac.uk/item/O221321/schotter-print-nees-georg/)
//! using whiskers.

use itertools::iproduct;
use whiskers::prelude::*;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::ops::{BitAnd, BitOr, Not};
use vsvg::{DocumentTrait, LayerTrait, DEFAULT_TOLERANCE};


/// Events for the scantable build. Each event is sorted by x, y, and slope in that order.
/// Ordering is done in *reversed* order to make the BinaryHeap structure give a minheap.
#[derive(Debug, Clone)]
pub struct Event {
    pub point: Point,
    pub index: i32,
    pub swap: Option<(i32, i32)>,
}

impl Eq for Event {}

impl PartialEq<Self> for Event {
    fn eq(&self, other: &Self) -> bool {
        return self.point.x == other.point.x && self.point.y == other.point.y;
    }
}

impl PartialOrd<Self> for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Event {
    fn cmp(&self, other: &Self) -> Ordering {
        match Point::cmp(&self.point, &other.point) {
            Ordering::Greater => {
                return Ordering::Less;
            }
            Ordering::Less => {
                return Ordering::Greater;
            }
            Ordering::Equal => {}
        }
        if other.index > self.index {
            return Ordering::Less;
        }
        if other.index < self.index {
            return Ordering::Greater;
        }
        return Ordering::Equal;
    }
}

#[derive(Debug, Clone)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

///Standard x,y point. Sort order for the point is x with tie breaks going to higher y-value.
impl Point {
    pub fn new(x: f64, y: f64) -> Point {
        Point { x, y }
    }
}

impl From<(f64, f64)> for Point {
    fn from(value: (f64, f64)) -> Self {
        Point::new(value.0, value.1)
    }
}

impl Eq for Point {}

impl PartialEq<Self> for Point {
    fn eq(&self, other: &Self) -> bool {
        return (self.x - other.x).abs() < 1e-12 && (self.y - other.y).abs() < 1e-12;
    }
}

impl PartialOrd<Self> for Point {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Point {
    fn cmp(&self, other: &Self) -> Ordering {
        if other == self {
            return Ordering::Equal;
        }
        if other.x < self.x {
            return Ordering::Greater;
        } else if other.x > self.x {
            return Ordering::Less;
        }
        if other.y < self.y {
            return Ordering::Greater;
        } else if other.y > self.y {
            return Ordering::Less;
        }
        return Ordering::Equal;
    }
}

/// Geomstr: Geometry class see, sister structure:
/// https://github.com/meerk40t/meerk40t/blob/main/meerk40t/tools/geomstr.py
#[derive(Debug, Clone)]
pub struct Geomstr {
    pub segments: Vec<((f64, f64), (f64, f64), (f64, f64), (f64, f64), (f64, f64))>,
}

impl Geomstr {
    pub fn new() -> Geomstr {
        Geomstr {
            segments: Vec::new(),
        }
    }
    pub fn from_segments(
        segments: Vec<((f64, f64), (f64, f64), (f64, f64), (f64, f64), (f64, f64))>,
    ) -> Geomstr {
        Geomstr { segments }
    }

    /// Add a rectangle to the geometry.
    pub fn rect(&mut self, x: f64, y: f64, width: f64, height: f64, settings: f64) {
        self.line((x, y), (x + width, y), settings);
        self.line((x + width, y), (x + width, y + height), settings);
        self.line((x + width, y + height), (x, y + height), settings);
        self.line((x, y + height), (x, y), settings);
    }

    /// Add a line to the geometry.
    pub fn line(&mut self, p0: (f64, f64), p1: (f64, f64), settings: f64) {
        self.segments
            .push((p0, (0., 0.), (41.0, settings), (0., 0.), p1));
    }

    /// Slope where divide by 0 is always negative infinity.
    pub fn slope(&self, index: usize) -> f64 {
        let line = &self.segments[index];
        let rise: f64 = line.0 .1 - line.4 .1;
        let run: f64 = line.0 .0 - line.4 .0;
        if run == 0.0 {
            return f64::NEG_INFINITY;
        }
        rise / run
    }

    /// Find an intersection between index0 and index1.
    pub fn get_intersection(&self, index0: usize, index1: usize) -> Option<(f64, f64)> {
        let line0 = &self.segments[index0];
        let line1 = &self.segments[index1];
        let a = &line0.0;
        let b = &line0.4;
        let c = &line1.0;
        let d = &line1.4;
        let denom: f64 = (d.1 - c.1) * (b.0 - a.0) - (d.0 - c.0) * (b.1 - a.1);
        if denom.abs() < 1e-12 {
            return None;
        }
        let t1: f64 = ((d.0 - c.0) * (a.1 - c.1) - (d.1 - c.1) * (a.0 - c.0)) / denom;
        let t2: f64 = ((b.0 - a.0) * (a.1 - c.1) - (b.1 - a.1) * (a.0 - c.0)) / denom;
        if 0.0 <= t1 && t1 <= 1.0 && 0.0 <= t2 && t2 <= 1.0 {
            return Some((t1, t2));
        }
        None
    }

    /// Returns the y_intercept point given a line a given x.
    /// Default is used for y if there is a line along the requested x.
    pub fn y_intercept(&self, index: usize, x: f64, default: f64) -> Point {
        let line = &self.segments[index];
        let a = &line.0;
        let b = &line.4;
        let rise: f64 = a.1 - b.1;
        let run: f64 = a.0 - b.0;
        if rise == 0.0 {
            return Point::new(x, a.1);
        }
        if run == 0.0 {
            return Point::new(x, default);
        }
        let m = run / rise;
        let x0: f64 = a.0 - (m * a.1);
        Point::new(x, (x - x0) / m)
    }

    /// Find point located within the current geometry at position t [0,1]
    pub fn point(&self, index: usize, t: f64) -> Point {
        let line = &self.segments[index];
        Point::new(
            t * (line.4 .0 - line.0 .0) + line.0 .0,
            t * (line.4 .1 - line.0 .1) + line.0 .1,
        )
    }
}

#[derive(Debug, Clone)]
pub struct BoolOp {
    pub inside: Vec<Vec<bool>>,
}

impl BoolOp {
    pub fn new(mask: Vec<Vec<bool>>) -> BoolOp {
        BoolOp { inside: mask }
    }
}

impl BitAnd for BoolOp {
    type Output = Self;

    fn bitand(self, rhs: Self) -> Self::Output {
        let mut n = Vec::new();
        for j in 0..rhs.inside.len() {
            let mut m = Vec::new();
            for k in 0..rhs.inside[j].len() {
                m.push(self.inside[j][k] & rhs.inside[j][k]);
            }
            n.push(m);
        }
        Self::new(n)
    }
}

impl BitOr for BoolOp {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        let mut n = Vec::new();
        for j in 0..rhs.inside.len() {
            let mut m = Vec::new();
            for k in 0..rhs.inside[j].len() {
                m.push(self.inside[j][k] | rhs.inside[j][k]);
            }
            n.push(m);
        }
        Self::new(n)
    }
}

impl Not for BoolOp {
    type Output = Self;

    fn not(self) -> Self::Output {
        let mut n = Vec::new();
        for j in 0..self.inside.len() {
            let mut m = Vec::new();
            for k in 0..self.inside[j].len() {
                m.push(!self.inside[j][k]);
            }
            n.push(m);
        }
        Self::new(n)
    }
}

#[derive(Debug, Clone)]
pub struct BeamTable {
    pub geometry: Geomstr,
    pub events: Vec<Point>,
    pub actives: Vec<Vec<i32>>,
    pub intersections: Vec<Point>,

    built: bool,
}

/// BeamTable acceleration structure. Creates a geometric space lookup table.
impl BeamTable {
    pub fn new(geometry: Geomstr) -> BeamTable {
        BeamTable {
            geometry,
            events: Vec::new(),
            actives: Vec::new(),
            intersections: Vec::new(),
            built: false,
        }
    }

    /// Create an Even/Odd fill for a given layer level.
    pub fn evenodd_fill(&self, layer: f64) -> BoolOp {
        let mut spacemask = Vec::new();
        for active in &self.actives {
            let mut active_mask = Vec::new();
            let mut inside = false;
            active_mask.push(inside);
            for a in active {
                let line = &self.geometry.segments[*a as usize];
                if line.2 .1 == layer {
                    inside = !inside;
                }
                active_mask.push(inside);
            }
            spacemask.push(active_mask);
        }
        BoolOp::new(spacemask)
    }

    /// Create an even_odd fill for all geometry.
    /// Useful for point in polygon solutions
    pub fn even_odd_ignoring_origin(&self) -> BoolOp {
        let mut spacemask = Vec::new();
        for active in &self.actives {
            let mut active_mask = Vec::new();
            let mut inside = false;
            active_mask.push(inside);
            for _a in active {
                inside = !inside;
                active_mask.push(inside);
            }
            spacemask.push(active_mask);
        }
        BoolOp::new(spacemask)
    }

    /// Create a union of all layers
    pub fn union_all(&self) -> BoolOp {
        let mut spacemask = Vec::new();
        for active in &self.actives {
            let mut set: HashMap<i32, bool> = HashMap::new();
            let mut active_mask = Vec::new();
            active_mask.push(set.len() != 0);
            for a in active {
                let line = &self.geometry.segments[*a as usize];
                if set.contains_key(&(line.2 .1 as i32)) {
                    set.remove(&(line.2 .1 as i32));
                } else {
                    set.insert(line.2 .1 as i32, true);
                }
                active_mask.push(set.len() != 0);
            }
            spacemask.push(active_mask);
        }
        BoolOp::new(spacemask)
    }

    /// Create geometry from a BoolOp.
    pub fn create(&self, mask: BoolOp) -> Geomstr {
        let mut g = Geomstr::new();
        let inside = &mask.inside;
        for j in 0..inside.len() - 1 {
            //mask exists at inside-1, but the final entry is actually pointless
            let left_event = &self.events[j];
            let beam_active = &self.actives[j];
            let right_event = &self.events[j + 1];

            for k in 0..inside[j].len() - 1 {
                let below_space = inside[j][k];
                let segment_active = beam_active[k];
                let above_space = inside[j][k + 1];
                if (below_space && !above_space) || (!below_space && above_space) {
                    //is a boundary.
                    let start =
                        self.geometry
                            .y_intercept(segment_active as usize, left_event.x, left_event.y);
                    let end =
                        self.geometry
                            .y_intercept(segment_active as usize, right_event.x, right_event.y);
                    let line = &self.geometry.segments[segment_active as usize];
                    g.line((start.x, start.y), (end.x, end.y), line.2 .1);
                }
            }
        }
        g
    }

    /// Find the actives for a particlar x/y event space.
    pub fn actives_at(&self, x: f64, y: f64) -> &Vec<i32> {
        let idx = self.events.binary_search(&Point::new(x, y));
        match idx {
            Ok(value) => {
                return &self.actives[value];
            }
            Err(value) => {
                if value == 0 {
                    return &self.actives.last().expect("at least 1 active must exist.");
                }
                let value = value.checked_sub(1).unwrap();
                return &self.actives[value];
            }
        }
    }

    /// Internal: find the position within the given actives for the current x.
    fn bisect_yints(&self, actives: &Vec<i32>, x: i32, scanline: &Point) -> i32 {
        let geometry = &self.geometry;
        let mut lo = 0;
        let mut hi = actives.len();
        let mut mid;
        while lo < hi {
            mid = (lo + hi) / 2;
            let test = &geometry.y_intercept(actives[mid] as usize, scanline.x, scanline.y);
            let value = &geometry.y_intercept(x as usize, scanline.x, scanline.y);
            match Point::cmp(&value, &test) {
                Ordering::Less => {
                    hi = mid;
                }
                Ordering::Greater => {
                    lo = mid + 1;
                }
                Ordering::Equal => {
                    let test_slope = &geometry.slope(actives[mid] as usize);
                    let value_slope = &geometry.slope(x as usize);
                    if value_slope < test_slope {
                        hi = mid
                    } else {
                        lo = mid + 1
                    }
                }
            }
        }
        lo as i32
    }

    /// Internal: check for intersections between indexes q and r, occurring after sl
    fn check_intersections(
        &mut self,
        events: &mut BinaryHeap<Event>,
        actives: &Vec<i32>,
        checked_swaps: &mut Vec<(i32, i32)>,
        q: usize,
        r: usize,
        sl: &Point,
    ) {
        let q = actives[q];
        let r = actives[r];
        let geometry = &self.geometry;
        if checked_swaps.contains(&(q, r)) {
            return;
        }
        let intersection = geometry.get_intersection(q as usize, r as usize);

        match intersection {
            None => (),
            Some(t) => {
                let t1 = t.0;
                let t2 = t.1;
                if (t1 == 0.0 || t1 == 1.0) && ((t2 == 0.0) || (t2 == 1.0)) {
                    return;
                }
                let pt_intersect = geometry.point(q as usize, t1);
                self.intersections.push(pt_intersect.clone());
                match Point::cmp(&sl, &pt_intersect) {
                    Ordering::Greater => {
                        return;
                    }
                    Ordering::Equal => {
                        return;
                    }
                    Ordering::Less => {}
                }
                checked_swaps.push((q, r));
                let event = Event {
                    point: pt_intersect,
                    index: 0,
                    swap: Some((q, r)),
                };
                events.push(event);
            }
        }
    }

    /// Builds the beamtable from the underlying geometry.
    pub fn build(&mut self) {
        if self.built {
            //This was already built.
            return;
        }
        let mut events: BinaryHeap<Event> = BinaryHeap::new();
        let mut checked_swaps: Vec<(i32, i32)> = Vec::new();
        let mut actives: Vec<i32> = Vec::new();

        // Create initial start and end values for the event queue.
        for i in 0..self.geometry.segments.len() {
            let line = &self.geometry.segments[i];
            let p0 = Point::new(line.0 .0, line.0 .1);
            let p1 = Point::new(line.4 .0, line.4 .1);
            match Point::cmp(&p0, &p1) {
                Ordering::Less => {
                    events.push(Event {
                        point: p0,
                        index: i as i32,
                        swap: None,
                    });
                    events.push(Event {
                        point: p1,
                        index: !i as i32,
                        swap: None,
                    });
                }
                _ => {
                    events.push(Event {
                        point: p1,
                        index: i as i32,
                        swap: None,
                    });
                    events.push(Event {
                        point: p0,
                        index: !i as i32,
                        swap: None,
                    });
                }
            }
        }

        // Process the event queue, performs Bentley-Ottmann line intersection checks
        while events.len() != 0 {
            let event = events
                .pop()
                .expect("Pop only called after checking events existed.");
            let idx = event.index;
            let index = event.index;
            let pt = &event.point;
            match event.swap {
                None => {
                    if idx >= 0 {
                        // Insert.
                        let ip = self.bisect_yints(&actives, index, &event.point) as usize;
                        actives.insert(ip, index);
                        if ip > 0 {
                            self.check_intersections(
                                &mut events,
                                &actives,
                                &mut checked_swaps,
                                ip - 1,
                                ip,
                                pt,
                            )
                        }
                        if ip < actives.len() - 1 {
                            self.check_intersections(
                                &mut events,
                                &actives,
                                &mut checked_swaps,
                                ip,
                                ip + 1,
                                pt,
                            )
                        }
                    } else {
                        //Remove.
                        let rp = actives
                            .iter()
                            .position(|&e| e == !index)
                            .expect("Was added should remove.");
                        actives.remove(rp);
                        if 0 < rp && rp < actives.len() {
                            self.check_intersections(
                                &mut events,
                                &actives,
                                &mut checked_swaps,
                                rp - 1,
                                rp,
                                pt,
                            )
                        }
                    }
                }
                Some((s1, _)) => {
                    let s1 = actives
                        .iter()
                        .position(|&e| e == s1)
                        .expect("Swap pos should exist.");
                    let s2 = s1 + 1;
                    actives.swap(s1, s2);
                    if s1 > 0 {
                        self.check_intersections(
                            &mut events,
                            &actives,
                            &mut checked_swaps,
                            s1 - 1,
                            s1,
                            pt,
                        );
                    }
                    if s2 < actives.len() - 1 {
                        self.check_intersections(
                            &mut events,
                            &actives,
                            &mut checked_swaps,
                            s2,
                            s2 + 1,
                            pt,
                        );
                    }
                }
            }
            match events.peek() {
                None => {}
                Some(last_pt) => {
                    if pt == &last_pt.point {
                        continue;
                    }
                }
            }

            // Push the current state to the table
            self.events.push((*pt).clone());
            self.actives.push(actives.clone());
        }
        self.built = true;
    }
}

#[sketch_app]
pub struct WhiskersDemoSketch {
    col_count: u32,
    row_count: u32,

    #[param(slider, min = 0., max = 10.)]
    offset_cm: f64,

    #[param(slider, min = 0., max = 10.)]
    box_size_cm: f64,

    #[param(slider, min = 0., max = 90.)]
    rand_angle_deg: f64,

    #[param(slider, min = 0., max = 3.)]
    rand_offset_cm: f64,

    #[param(slider, min = 0., max = 10.)]
    stroke_width: f64,
}

impl Default for WhiskersDemoSketch {
    fn default() -> Self {
        Self {
            col_count: 12,
            row_count: 24,
            offset_cm: 1.,
            box_size_cm: 1.,
            rand_angle_deg: 45.,
            rand_offset_cm: 0.3,
            stroke_width: 1.0,
        }
    }
}

impl App for WhiskersDemoSketch {
    fn update(&mut self, sketch: &mut Sketch, ctx: &mut Context) -> anyhow::Result<()> {
        sketch.scale(Unit::Cm).stroke_width(self.stroke_width);

        for (i, j) in iproduct!(0..self.col_count, 0..self.row_count) {
            sketch.push_matrix_and(|sketch| {
                sketch.translate(i as f64 * self.offset_cm, j as f64 * self.offset_cm);

                let max_angle = self.rand_angle_deg * (j as f64 / self.row_count as f64);
                let max_offset = self.rand_offset_cm * (j as f64 / self.row_count as f64);

                sketch
                    .rotate_deg(ctx.rng_range(-max_angle..max_angle))
                    .translate(
                        ctx.rng_range(-max_offset..max_offset),
                        ctx.rng_range(-max_offset..max_offset),
                    )
                    .rect(0., 0., self.box_size_cm, self.box_size_cm);
            });
        }
        let doc = sketch.document_mut();
        let doc = doc.flatten(DEFAULT_TOLERANCE);

        // convert everything to lines
        let mut segments = Geomstr::new();
        let mut idx = 0;
        doc.layers.values().for_each(|layer| {
            layer.paths.iter().for_each(|path| {
                path.data.points().windows(2).for_each(|p| {
                    segments.line((p[0].x(), p[0].y()), (p[1].x(), p[1].y()), idx as f64);
                });
                idx += 1;
            });
        });

        // run scan beam algorithm
        let mut beamtable = BeamTable::new(segments);
        beamtable.build();
        // let mask = beamtable.evenodd_fill(20.0);
        let mask = beamtable.union_all();
        let geom = beamtable.create(mask);

        //
        // visualize the result
        //

        // convert back to regular (not flattened) document, merge everything to layer 0 and normalize
        // line width and color
        // let mut new_doc = vsvg::Document::default();

        let layer = doc.get_mut(1);
        for line in geom.segments {
            layer.line(line.0 .0, line.0 .1, line.4 .0, line.4 .1);
        }

        doc.merge_layers();
        doc.for_each(|layer| {
            layer.for_each(|path| {
                path.metadata_mut().stroke_width = 0.5;
                path.metadata_mut().color = vsvg::Color::LIGHT_RED;
            });
        });
        // *doc = new_doc;

        Ok(())
    }
}

wasm_sketch!(WhiskersDemoSketch::runner()
    .with_layout_options(LayoutOptions::centered())
    .with_info_options(
        InfoOptions::default()
            .description(
                "This sketch is a recreation of the classic \"Schotter\" series by Georg Nees \
            (1968-1970).\n\nGeorg Nees (born 1926, Nuremberg) is considered one of the founders \
            of computer art and graphics. He was also one of the first people to exhibit his \
            computer graphics, at the studio gallery of the Technische Hochschule in Stuttgart in \
            February 1965. In 1969, he received his doctorate on the subject of Generative \
            Computer Graphics."
            )
            .author("Antoine Beyeler")
            .author_url("https://bylr.info/")
            .source_url(
                "https://github.com/abey79/vsvg/blob/master/crates/whiskers-web-demo/src/lib.rs"
            )
    ));
