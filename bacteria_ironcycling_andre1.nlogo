;; first bacteria model (universal model for further ChemBioSys projects)

;;;;;;;;;;;;;;;;;
;; inits of individuals and variables

breed [ fe3reducer fe3reducers ]
; breed [ fe3reducerB fe3reducersB ]  ;; possibility to split into fe3reducer (A) and fe3reducerB
breed [ fe2oxidizer fe2oxidizers ]


turtles-own [
  energy ;; fitness of each individual
  ;partner ;; interaction parter (to exchange products/information/...)
  ;substrate-essential  ;; substrate/nutrient that is primarily used/processed (e.g. an essential amino acid, or a main nutrient); taken from other breed
  ;substrate-dispensable  ;; substrate/nutrient that is produced or taken up from the environment; can be transfered to other breed
  split-ticks ;; tick of reproduction (split into two bacteria)
  in-biofilm? ;; True if they are part of a biofilm, False if they are not
]

patches-own [
  fe?   ;; Boolean, True if the field contains iron
  redox-balance   ;; negative values if there has been more reduction on this field (Fe2+), positive values if more oxidation (Fe3+)
  cluster-identifier   ;; float identifier shared by iron containing patches of the same cluster, 1 if they don`t contain iron
]

globals[
  max-ticks  ;; simulation runs for one hour (6 ticks = 1 minute)
  default-color

  ;; single descriptive values for plotting
  global-mean-energy-fe3reducer
  ; global-mean-energy-fe3reducerB
  global-mean-energy-fe2oxidizer

  global-reduction-count  ;; count the number of reduction reactions
  global-oxidation-count  ;; count the number of oxidation reactions
  global-redox-count ;; global redox balance

  ;; lists and variables needed for the move-shewanella function
  cluster-identifier-list ; list of identifiers of clusters close to fe3-reducer
  cluster-redox-balance-list ; list of arithmetic means of redox-balances of each patch for each cluster
  neighborcluster-total-redox-balance ; arithmetic mean of cluster-redox-balance-list's elements (overall redox-balance of patches nearby)
]

;;;;;;;;;;;;;;;;;
;; initial setup of simulation

to setup
  clear-all

  ;; color patches
  ask patches [
    set pcolor white
    set cluster-identifier 1
  ]

  draw-walls  ;; create borders on all sites

  ;; create fe-containing patches randomly
  ;; non-clustered fe-patches
  repeat ((count patches - count patches with [pcolor = black]) * (fe-patch-percentage / 100) * (1 - fe-clustering / 100)) [
    ifelse count patches with [(pcolor != black) and
                  (fe? != True) and
                  (not (any? neighbors with [ fe? = True ]))] > 0 [
      assign-fe False
    ]
    [
      assign-fe True
    ]
  ]
  ;; clustered fe-patches
  repeat ((count patches - count patches with [pcolor = black]) * (fe-patch-percentage / 100) * (fe-clustering / 100)) [
    ifelse count patches with [ fe? = True ] > 0 [
      assign-fe True
    ]
    [
      assign-fe False
    ]
  ]

  update-color
  while [ count patches with [(fe? = True) and (cluster-identifier = 1)] > 0]
    [ask n-of 1 (patches with [(fe? = True) and (cluster-identifier = 1)])
      [ set cluster-identifier random-float 1
        group-clusters ]
    ]

  set-default-shape turtles "checker piece"

  set max-ticks 3700 ;; run it over a longer time period?

  ;; create bacteria species
  create-fe3reducer initial-number-fe3reducer [
    set color blue
    set size 1
    set in-biofilm? False
    set energy 1 + random fe3reducer-max-initial-energy  ;; energy is not normal distributed because all individuals start at a different state, not all with full energy
    randomize

    ;; define ticks of each species to reproduce (same distance for all species but different time point)
    let split-time-range max-ticks / doubling-X-times-fe3reducer  ;; number of cell divisions of each bacterium (see Schaechter1962, Table 1 for some values)
    let split-tick1 random split-time-range
    set split-ticks []
    let i 0
    while [i < max-ticks and split-tick1 < max-ticks][
      set split-ticks lput (split-tick1 + i * split-time-range) split-ticks
      set i i + 1
    ]
  ]

  create-fe2oxidizer initial-number-fe2oxidizer [
    set color red
    set size 1
    set in-biofilm? False
    set energy 1 + random fe2oxidizer-max-initial-energy  ;; energy is not normal distributed because all individuals start at a different state, not all with full energy
    randomize

    ;; define ticks of each species to reproduce (same distance for all species but different time point)
    let split-time-range max-ticks / doubling-X-times-fe2oxidizer
    let split-tick1 random split-time-range
    set split-ticks []
    let i 0
    while [i < max-ticks and split-tick1 < max-ticks][
      set split-ticks lput (split-tick1 + i * split-time-range) split-ticks
      set i i + 1
    ]
  ]

  reset-ticks
  update-reporter
  set global-reduction-count 0
  set global-oxidation-count 0
end


;;;;;;;;;;;;;;;;;
;; runtime procedures, for each time step

to go
  if not any? turtles [ stop ]
  if ticks = max-ticks + 1 [ stop ]
  switch-oxygen

  ask turtles [
    if show-labels? [ set label "" ]
    if member? ticks split-ticks and energy > 5 [  ;; minimum energy necessary for reproduction
      reproduce  ;reproduce - depending on doubling-X-times (defined in split-ticks)
    ]
  ]

  if O2-saturation > 10 [
    fe2-oxidation-by-air
  ]

  ask fe3reducer [ ;; move, die, reduce Fe3 to Fe2
    if (in-biofilm? != True)[
      move-shewanella speed-fe3reducer
    ]
    death

    ;; check concentration of Fe3 on patch under agent
    ifelse O2-saturation < 5 [
      if [ fe? = True ] of patch-here [
        ;; Shewanella uses O2 as an electron sink
        ;; only when O2 is missing, it reduces Fe3 to Fe2 to gain energy
        if (random-float 1) > (-1 * [redox-balance] of patch-here * passivation-factor) [
        ;; only when there is no Fe2+ on the surface, Fe3+ can be reduced to Fe2+
          ask patch-here [
            set redox-balance redox-balance - 1
          ]
          set energy energy + 1
          if show-labels? [ set label "R" ]
          set global-reduction-count global-reduction-count + 1
        ]
      ]
    ]
    [
      ;; get energy from O2, approximately five times the amount of iron respiration (Kostka2002).
      set energy energy + 5
    ]
  ]

  ask fe2oxidizer [ ;; move, die, oxidize Fe2 to Fe3
    if (in-biofilm? != True)[
      move-sideroxydans speed-fe2oxidizer
    ]
    death

    ;; check concentration of Fe3 on patch under agent
    ifelse O2-saturation > 17 and  O2-saturation < 21 [
      if [redox-balance] of patch-here < 0 [
        ;; Sideroxydans uses O2 as energy source
        ;; only between 18-20% oxidizes Fe2 to Fe3
        ask patch-here [
          set redox-balance redox-balance + 1
        ]
        set energy energy + 5
        if show-labels? [ set label "O" ]
        set global-oxidation-count global-oxidation-count + 1
      ]
    ]
    [
      ;; get energy from alternative electron sink, e.g. nitrate
      set energy energy + 1
    ]
  ]

  make-end-biofilms
  update-color
  update-reporter

end



;;;;;;;;;;;;;;;;;
;; observer procedure that creates borders on all sites
to draw-walls
  ;; don't let the window be bigger than the right chamber
  ask patches with [(pxcor = min-pxcor) or
                    ((pxcor < 0) and (abs pycor = max-pycor)) or
                    (pxcor >= max-pxcor) or
                    ((pxcor > 0) and (abs pycor >= max-pycor))]
    [ set pcolor black ]
  ask patches with [pxcor = 0] [
    ifelse abs pycor < 20  ;; if it is smaller than 20, there is a wall between right and left side (see Simple Kinetics 3)
        [ ] ;set pcolor brown + 1]
        [ set pcolor black ]
  ]

  ;; if the window size changed
  ask turtles with [(pxcor = 0) and (pcolor = black)]
    [ randomize ]
end

;;switches oxygen-levels to specified points of time.
to switch-oxygen
  let oxy-on (range 600 ((600 * 6) + 1) 600)
  let oxy-off (range 700 ((700 * 6) + 1) 600)
  if (member? ticks oxy-on) = True [
    set O2-saturation 20
  ]
  if (member? ticks oxy-off) = True [
    set O2-saturation 3
  ]
end

;;;;;;;;;;;;;;;;;
;; patch procedures
to assign-fe [ cluster? ]
  ifelse cluster? = False [
    ask n-of 1 (patches with [(pcolor != black) and
                  (fe? != True) and
                  (not (any? neighbors with [fe? = True]))]) [
      set pcolor red
      set fe? True
      set redox-balance 0
    ]
  ]
  [ ask n-of 1 (patches with [(pcolor != black) and
               (fe? = True) and
               (count neighbors with [(fe? = True) or (pcolor = black)] < 8)]) [
        ask (last sort-on [count neighbors with [fe? = True]] neighbors with [(fe? != True) and (pcolor != black)]) [
          set pcolor red
          set fe? True
          set redox-balance 0
        ]
    ]
  ]
end

to update-color
  ask patches with [(pxcor != min-pxcor) and
                    ((pxcor > 0) or (abs pycor != max-pycor)) and
                    (pxcor != max-pxcor) and
                    ((pxcor <= 0) or (abs pycor < max-pycor)) and
                    Fe? = True][
    ifelse redox-balance = 0
     [set pcolor red]  ;; neutral balance
     [
      ifelse redox-balance < 0
       [set pcolor blue + 2 + (redox-balance * passivation-factor)]  ;; more Fe2+ produced recently
       [set pcolor red + 2 - (redox-balance * passivation-factor)]  ;; more Fe3+ produced recently
     ]
  ]
end

;;assigns same cluster-identifier to patches of the same cluster
to group-clusters
  ask neighbors4 with [(cluster-identifier = 1) and
    (fe? = True)]
  [ set cluster-identifier [cluster-identifier] of myself
    group-clusters ]
end

;;fe2 is oxidized through oxygen
to fe2-oxidation-by-air
  ask patches with [fe? = True] [
    if redox-balance < 0 and random-float 1 > 0.5 [
      set redox-balance redox-balance + 1
      set global-oxidation-count global-oxidation-count + 1
    ]
  ]
end
;;;;;;;;;;;;;;;;;
;; turtle procedures

to randomize
  setxy random-xcor random-ycor
  if any? patches in-radius 1 with [pcolor = black]
    [ randomize ] ;; keep trying until we don't land on or near black
end

to move-shewanella [ speed ]  ;; move turtles around randomly ;; let them stop if a wall is ahead
  ;; if patch-here has fe the turtle will reverse with probability 1, if not it will reverse with probability of 1/3
  ifelse ((any? ((patch-set patch-here neighbors) with [fe? = True])) and 
          (mean [redox-balance] of ((patch-set patch-here neighbors) with [fe? = True]) >= (-0.5 / passivation-factor))
  [
    rt  random 360
    ifelse [pcolor] of patch-ahead 1 != black
    [
      fd speed * 1
    ]
    [
      move-shewanella speed
    ]
  ]
  [
    if random-float 1 <= (1 / 3) [rt random 360]
    ifelse [pcolor] of patch-ahead 1 != black
    [
      fd speed * 0.5
    ]
    [
      move-shewanella speed
    ]
  ]

to move-sideroxydans [ speed ]  ;; move turtles around randomly ;; let them stop if a wall is ahead
  rt random 360
  ifelse [pcolor] of patch-ahead 1 != black
  [
    fd speed * 0.5 ;; TODO adjust
  ]
  [
    move-shewanella speed
  ]
end


to death
  ;; die if you run out of energy
  if energy <= 0 [ die ]
end


to reproduce
  ;; devide into two bacteria, split energy to both individuals
  let new-energy energy / 2

  ;; define new for parent
  set energy new-energy

  let children-in-biofilm? in-biofilm?

  hatch 1 [
    set energy new-energy
    set split-ticks n-values 3 [random 0]  ; do not reproduce for now
    set in-biofilm? children-in-biofilm?
  ]

end

;; create, expand or resolve biofilms.
to make-end-biofilms
  ;; only adress iron-cluster patches.
  let unique-clusters remove-duplicates [cluster-identifier] of patches with [cluster-identifier < 1]
  ;; adress every cluster.
  foreach unique-clusters [ x ->
    ;; There is no biofilm on this cluster yet.
    ifelse (not (any? (fe3reducer-on (patches with [cluster-identifier = x])) with [in-biofilm? = True]))[
      ;; check if the cluster is big enough, there are enough bacteria on the patch and suffers not too much passivation.
      if ((count patches with [cluster-identifier = x] >= 4) and
        (count fe3reducer-on patches with [cluster-identifier = x] >= 4) and
        ((mean [redox-balance] of patches with [cluster-identifier = x]) * passivation-factor * -1 <= biofilm-start-treshold)) [
        ask fe3reducer-on patches with [cluster-identifier = x] [
          set in-biofilm? True
          let patches-of-cluster patches with [cluster-identifier = x]
          move-to (max-one-of (patches-of-cluster with-min [count fe3reducer-here]) [count neighbors with [any? fe3reducer-here with [in-biofilm? = True]]])
        ]
      ]
    ]
    ;; there is already a biofilm on the cluster.
    [
      ;; if the passivation is too high, resolve the biofilm stepwise.
      ifelse ((mean [redox-balance] of patches with [cluster-identifier = x]) * passivation-factor * -1 >= biofilm-dissolve-treshold) [
        ;; if the biofilm is already quite small, release all the bacteria.
        ifelse (count fe3reducer-on patches with [cluster-identifier = x] <= 4)[
          ask fe3reducer-on patches with [cluster-identifier = x] [
            set in-biofilm? False
          ]
        ]
        ;; release only some of the bacteria
        [
          ask fe3reducer-on patches with [cluster-identifier = x] [
            if (random-float 1 <= biofilm-release-probability) [
              set in-biofilm? False
            ]
          ]
        ]
      ]
      ;; if the passivation is still quite low, allow recruitment of new bacteria
      [
        ask (fe3reducer-on (patches with [cluster-identifier = x])) with [in-biofilm? = False][
          if (random-float 1 <= biofilm-recruitment-probability)[
            set in-biofilm? True
            let patches-of-cluster patches with [cluster-identifier = x]
            move-to (max-one-of (patches-of-cluster with-min [count fe3reducer-here]) [count neighbors with [any? fe3reducer-here with [in-biofilm? = True]]])
          ]
        ]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;
;; update plots

to update-reporter
  ;; update global counter for plots
  ifelse any? fe3reducer
  [ set global-mean-energy-fe3reducer mean [ energy ] of fe3reducer ]
  [ set global-mean-energy-fe3reducer 0 ]

;  ifelse any? fe3reducerB
;  [ set global-mean-energy-fe3reducerB mean [ energy ] of fe3reducerB ]
;  [ set global-mean-energy-fe3reducerB 0 ]


  ifelse any? fe2oxidizer
  [ set global-mean-energy-fe2oxidizer mean [ energy ] of fe2oxidizer ]
  [ set global-mean-energy-fe2oxidizer 0]


  ifelse any? turtles and ticks > 0
  [
    set global-redox-count global-reduction-count - global-oxidation-count
  ]
  [
    set global-redox-count 0
  ]

end

;;;;;;;;;;;;;;;;;
@#$#@#$#@
GRAPHICS-WINDOW
361
66
900
606
-1
-1
12.95122
1
10
1
1
1
0
1
1
1
-20
20
-20
20
1
1
1
ticks
30.0

BUTTON
18
19
91
52
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
21
178
287
211
fe3reducer-max-initial-energy
fe3reducer-max-initial-energy
0
100
47.0
1
1
NIL
HORIZONTAL

BUTTON
109
18
172
51
go
go\ntick
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
24
450
287
483
fe2oxidizer-max-initial-energy
fe2oxidizer-max-initial-energy
0
100
62.0
1
1
NIL
HORIZONTAL

PLOT
920
80
1370
256
Energy
time
energy-level
0.0
400.0
0.0
50.0
true
true
"" ""
PENS
"energy fe3reducer" 1.0 0 -7500403 true "" "if count fe3reducer > 0 [plot global-mean-energy-fe3reducer]\n"
"energy fe2oxidizer" 1.0 0 -2674135 true "" "if count fe2oxidizer > 0 [plot global-mean-energy-fe2oxidizer]"

SLIDER
21
217
211
250
speed-fe3reducer
speed-fe3reducer
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
23
490
210
523
speed-fe2oxidizer
speed-fe2oxidizer
0
1
0.5
0.01
1
NIL
HORIZONTAL

PLOT
921
262
1372
412
Count
time
count
0.0
400.0
0.0
60.0
true
true
"" ""
PENS
"number fe3reducer" 1.0 0 -7500403 true "" "plot count fe3reducer"
"number fe2oxidizer" 1.0 0 -2674135 true "" "plot count fe2oxidizer"

INPUTBOX
21
113
182
173
initial-number-fe3reducer
200.0
1
0
Number

INPUTBOX
23
384
184
444
initial-number-fe2oxidizer
200.0
1
0
Number

TEXTBOX
24
343
229
383
setup Fe(II) oxidizing bacteria;\nhere Sideroxydans CL21\n
12
15.0
1

TEXTBOX
22
73
238
112
setup Fe(III) reducing bacteria; \nhere Shewanella oneidensis MR-1
12
5.0
1

TEXTBOX
927
17
1188
66
1 tick represents 10 seconds,\n6 ticks represent 1 minute,\n360 ticks represent 1 hour.
13
0.0
1

SWITCH
697
652
845
685
show-labels?
show-labels?
0
1
-1000

TEXTBOX
368
20
656
50
Plate size represents 160 µm = 0.16 mm\none bacterium has a size of 2 µm (default)
12
0.0
1

PLOT
920
417
1373
567
Redox-Balance
time
Fe2+
0.0
400.0
0.0
1000.0
true
true
"" ""
PENS
"Fe-redox-balance" 1.0 0 -16777216 true "" "plot global-redox-count"

TEXTBOX
930
580
1143
625
patch color: concentration of Fe(II) or/and Fe(III) > 0
12
0.0
1

TEXTBOX
932
629
1082
647
Fe(III)
12
17.0
1

TEXTBOX
933
652
1083
670
Fe(II) on surface
12
107.0
1

TEXTBOX
933
672
1083
690
otherwise
12
36.0
1

TEXTBOX
1392
632
1600
670
setup Fe(III) reducing bacteria;\nhere Pseudomonas fluorescens\n
12
55.0
1

INPUTBOX
1390
668
1551
728
initial-number-fe3reducerB
100.0
1
0
Number

SLIDER
1389
735
1665
768
fe3reducerB-max-initial-energy
fe3reducerB-max-initial-energy
0
50
35.0
1
1
NIL
HORIZONTAL

SLIDER
1389
773
1589
806
speed-fe3reducerB
speed-fe3reducerB
0
1
0.2
0.01
1
NIL
HORIZONTAL

SLIDER
1389
812
1724
845
doubling-X-times-fe3reducerB
doubling-X-times-fe3reducerB
0
2
0.1
0.1
1
per h
HORIZONTAL

SLIDER
361
614
615
647
fe-patch-percentage
fe-patch-percentage
0
100
20.0
1
1
%
HORIZONTAL

INPUTBOX
23
532
184
592
doubling-X-times-fe2oxidizer
0.694
1
0
Number

INPUTBOX
20
256
181
316
doubling-X-times-fe3reducer
1.453
1
0
Number

SLIDER
361
657
615
690
fe-clustering
fe-clustering
0
100
92.0
1
1
%
HORIZONTAL

SLIDER
362
701
615
734
passivation-factor
passivation-factor
0
1
0.02
0.02
1
NIL
HORIZONTAL

SLIDER
697
612
870
645
O2-saturation
O2-saturation
0
100
3.0
1
1
%
HORIZONTAL

SLIDER
25
636
197
669
biofilm-start-treshold
biofilm-start-treshold
0
1
0.2
0.05
1
NIL
HORIZONTAL

SLIDER
25
678
212
711
biofilm-dissolve-treshold
biofilm-dissolve-treshold
0
1
0.5
0.05
1
NIL
HORIZONTAL

SLIDER
25
720
221
753
biofilm-release-probability
biofilm-release-probability
0
1
0.25
0.05
1
NIL
HORIZONTAL

SLIDER
25
762
244
795
biofilm-recruitment-probability
biofilm-recruitment-probability
0
1
0.5
0.05
1
NIL
HORIZONTAL

@#$#@#$#@
## Purpose

The models purpose is to gain insights from interplay of iron reducing (e.g. Shewanella oneidensis MR-1) and iron oxidizing bacteria (e.g. Sideroxydans CR21). The goal is to predict the optimal O2-dynamics for both bacteria in terms of energetics/reproduction..

## Entities, state variables, and scales

We start off with two agent breeds, one iron reducing and one iron oxidizing bacteria species. The patches either contain iron or not. The amount of iron patches and how clustered they are can be definied via the sliders. As the bacteria reduce the iron the initially red patch starts of light blue and gets more darker while the reduced iron (Fe2+) accumulates. 6 ticks represent 1 minute and therefore 360 ticks equal 1 hour. One bacteria is 2 µm in diameter and the plate has an edge length of 0.16 mm. For the bacteria their number, maximum initial energy, speed and doubling time can be defined. The passivation-factor describes how much harder it becomes for Shewanella to reduce iron once the surface is covered with more and more Fe2+ (High passivation factor means that it gets easier to reduce iron). Additionally O2-saturation of the environment may be defined.

## Process overview and scheduling



## Design concepts



## Initialization



## Input data



## Submodels



## Credits and references

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

checker piece
false
0
Circle -7500403 true true 60 60 180
Circle -16777216 false false 60 60 180
Circle -7500403 true true 75 45 180
Circle -16777216 false false 75 45 180

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

circlebac
false
4
Circle -7500403 true false 0 0 300

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

egg
false
0
Circle -7500403 true true 96 76 108
Circle -7500403 true true 72 104 156
Polygon -7500403 true true 221 149 195 101 106 99 80 148

ellipse
true
2
Circle -7500403 true false 30 90 88
Circle -7500403 true false 45 75 90
Circle -7500403 true false 60 60 90

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

spinner
true
0
Polygon -7500403 true true 150 0 105 75 195 75
Polygon -7500403 true true 135 74 135 150 139 159 147 164 154 164 161 159 165 151 165 74

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="exampleRun" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="360"/>
    <metric>count fe3reducer</metric>
    <metric>count fe3reducerB</metric>
    <metric>count fe2oxidizer</metric>
    <metric>global-mean-energy-fe3reducer</metric>
    <metric>global-mean-energy-fe3reducerB</metric>
    <metric>global-mean-energy-fe2oxidizer</metric>
    <metric>global-oxidation-rate</metric>
    <metric>global-reduction-rate</metric>
    <enumeratedValueSet variable="doubling-X-times-fe2oxidizer">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fe2oxidizer-max-initial-energy">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-labels?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overall-fe-concentration">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number-fe3reducer">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number-fe3reducerB">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="speed-fe3reducerB">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number-fe2oxidizer">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="speed-fe2oxidizer">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion-fe3-to-fe2">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="speed-fe3reducer">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fe3reducerB-max-initial-energy">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doubling-X-times-fe3reducerB">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fe3reducer-max-initial-energy">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doubling-X-times-fe3reducer">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
