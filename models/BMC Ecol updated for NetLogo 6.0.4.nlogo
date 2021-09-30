;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Netlogo model of competition and colonisation in a spatially explicit environment
; A. L. Cronin, N. Loeuille & T. Monnin
; March 2015
;
; This model investigates the success of different investment strategies (competition
; /colonization tradeoff) under varying environmental conditions, with a specific focus
; on Independent (better coloniser) and Dependent (better competitor) colony foundation in
; social insects.
;
; All mature agents (referred to as turtles in NetLogo) can produce offspring. This assumes that
; they are all females (i.e. a queen heading a colony) and that males are not a limiting factor
; and only provide mating. In particular, males provide no help in collecting resources and
; producing offspring



;;;;;;;;;; DEFINITION OF VARIABLES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

breed [s1-agents s1-agent]   ; defines two categories of agents (strategy 1 and strategy 2)
breed [s2-agents s2-agent]


patches-own
[
  pclass                     ; patches belong to one of two classes, class 1 and class 2 patches
  pquality                   ; long term quality of the patch = mean of presources at each tick
  presources                 ; short term quality of the patch = actual presources available
]


turtles-own
[
  tresources                 ; amount of resources agents collect from their patch and use for
                             ; maintenance, growth and reproduction (i.e. body size / worker
                             ; number + stored energy)
  age                        ; age of agents
  death-age                  ; age at which agents die of old age
  maintenance-cost           ; cost incurred by agents to maintain themselves
  max-growth-rate            ; determines the growth rate of agents when they are not limited by
                             ; resources
  maturity-threshold         ; resources above which agents reproduce (i.e. sexual maturity)
  reproductive-investment    ; percent of tresources invested in reproduction

  offspring-number           ; number of offspring produced
  dispersal-range            ; maximal distance at which a new agent may disperse from her mother
  dispersal-risk             ; risk of mortality during dispersal (= % of offspring that die while
                             ; dispersing)
  dispersal-patch-choice     ; criteria for selecting dispersal site
]



;;;;;;;;;; SETUP ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all
  setup-patches
  update-displayed-values
  setup-turtles
  reset-ticks
end



;;;;;;;;;; GO ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go
  update-presources               ; replenishes resources in patches according to the
                                  ; characteristics of each patch
  update-displayed-values
    ; updating displays here allows checking that presources are correctly RESET at each tick,
    ; including that variance is correctly implemented

  partition-presources
    ; agents sharing a patch partition presources. Partitioning depends on relative tresources
    ; (i.e. body size / worker number + stored energy)

  maintain-turtle                 ; agents consume tresources to maintain themselves. If they
                                  ; don't acquire enough resources they shrink
  reproduce                       ; agents with more than a threshold resource level reproduce
  live-or-die                     ; agents may die from insufficient resources (starvation),
                                  ; old age, or random events

  if invasion? and ticks = invasion-time [launch-invasion]
    ; used to study invasion of a resident strategy at equilibrium by the other strategy

  if stop-on-win? and check-any-strategy-dead = true [stop]
    ; stop the simulation when one strategy has been eradicated by the other

  tick
end


;;;;;;;;;; SETUP LANDSCAPE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-patches
; Generates the landscape of patches. Generates three types of landscape with varied pclass
; and pquality :
; - a uniform landscape composed of patches of medium quality
; - a patchy landscape composed of 50% good and 50% bad patches that are randomly distributed
; - a patchy landscape composed of 50% good and 50% bad patches that are  more or less aggregated
; (3 levels of aggregation)
; The fractal component used for generating the aggregated landscapes is based on a fractal
; generator by James Steiner

    ifelse landscape-aggregation = "uniform"
    [                                ; creates a continuous 'uniform' landscape with pquality
       ; scaled to the mean of settings for patchy landscapes, so that overall mean resource
       ; input is the same
      ask patches
      [
       ifelse random-float 1 < prop-class1-patches
       [pclass1]                     ; each patch has a probability of being of class 1 equal to
                                     ; slider-defined landscape-wide proportion of class 1 patches
       [pclass2]                     ; patches not class 1 are class 2
       set pquality (prop-class1-patches + ((1 - prop-class1-patches) *
           class2-patches-relative-pquality)) ; make pquality uniform across patches
       set pcolor 57                 ; intermediate green to signify difference
                                     ; from full resource input
      ]

    ]
    [                                ; if the landscape is not uniform check whether
                                     ; it is random
      ifelse landscape-aggregation = "random"
      [
        ask patches
        [
          ifelse random-float 1 < prop-class1-patches
          [pclass1]
          [pclass2]
        ]
      ]
      [       ; if the landscape is neither uniform nor random
              ; then it's fractal. Generates a landscape of more or less aggregated patches
        let variation 10
        let dim min (list world-width world-height)
        let g 2 ^ ( floor log dim 2 - 1 )
        let diamond-set 0
        let square-set 0
        ask patches [set pquality 0]
        set diamond-set (patch-set patch 0 0)
        repeat floor log world-width 2
        [
          set square-set (patch-set [my-parent-corners g] of diamond-set)
          ask square-set [set pquality (pquality + (1 - random-float 2) * variation
              + mean [pquality] of my-parent-corners g)]
          set diamond-set (patch-set [my-parent-points g] of square-set)
          ask diamond-set [set pquality (pquality + (1 - random-float 2) * variation
              + mean [pquality] of my-parent-points g)]
          set g g / 2
          set variation variation * (2 ^ (- landscape-aggregation))
        ]

        let HQ-patch-set max-n-of (int (count patches * prop-class1-patches)) patches [pquality]
          ; retains the slider-defined number of patches with highest pquality
        ask patches
        [
          ifelse member? self HQ-patch-set
          [pclass1]                    ; patches that belong to this set of patches are
                                       ; class 1 patches
          [pclass2]                    ; the remaining patches are class 2
        ]
      ]
    ]

    update-presources                  ; set up is completed by generating presources from
                                       ; pquality, including implementing ± 10% variance
end



;;;;;;;;;; REPORTS & SUB-PROCEDURES OF SETUP LANDSCAPE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report my-parent-corners [g]
  report patches at-points (list (list (- g) g) (list g g) (list g (- g)) (list (- g) (- g) ) )
end


to-report my-parent-points [g]
  report patches at-points (list (list 0 (- g) ) (list g 0) (list 0 g) (list (- g) 0) )
end


to pclass1
; class 1 patches are the good patches. They have pquality = 1 by definition and are shown
; in green
  set pclass 1
  set pquality 1
  set pcolor green
end


to pclass2
; class 2 patches are the bad patches. They have slider-defined pquality (between 0 and 1) and
; are shown in brown
  set pclass 2
  set pquality class2-patches-relative-pquality
  set pcolor brown
end


; Note that, in the uniform landscape, patches of medium quality are not coded as such (see ifelse
; landscape-aggregation = "uniform" above).


to update-displayed-values
; this allows visual checks of pquality and presources spatial distribution and temporal variation

  ask patches
  [
    if show-pquality-value? = false and show-presources-value? = false [set plabel ""]
    if show-pquality-value? [set plabel precision pquality 2]
    if show-presources-value? [set plabel presources]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;; SETUP AGENTS

to setup-turtles
  set-default-shape turtles "circle"

  create-s1-agents start-agents-s1
  [
    set tresources (1 + random (start-resources - 1))
      ; sets tresources between 1 and start-resources
    set age random longevity
    set death-age longevity
    set maintenance-cost maintenance
    set max-growth-rate max-growth
    set maturity-threshold maturity-age
    set reproductive-investment investment-in-reproduction

    set offspring-number offspring-number-s1
    set dispersal-range dispersal-range-s1
    set dispersal-risk dispersal-risk-s1
    set dispersal-patch-choice dispersal-patch-choice-s1

    setxy random-xcor random-ycor
    set color red
  ]

  create-s2-agents start-agents-s2
  [
    set tresources (1 + random (start-resources - 1))
    set age random longevity
    set death-age longevity
    set maintenance-cost maintenance
    set max-growth-rate max-growth
    set maturity-threshold maturity-age
    set reproductive-investment investment-in-reproduction

    set offspring-number offspring-number-s2
    set dispersal-range dispersal-range-s2
    set dispersal-risk dispersal-risk-s2
    set dispersal-patch-choice dispersal-patch-choice-s2

    setxy random-xcor random-ycor
    set color blue
  ]

  ifelse show-turtle-size?
    ; when selected, allows visualising distribution of agent size across the landscape,
    ; and changes in agent size over time
  [update-turtle-size]
  [ask turtles [set size 0.3]]
end



;;;;;;;;;; SUB-PROCEDURES OF SETUP AGENTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-turtle-size
  ask turtles
  [
    set size (tresources / maturity-age)   ; scales size of turtles icons proportional to
                                           ; tresources
    if size < 0.15 [set size 0.15]         ; very small turtles are shown larger than
                                           ; proportionally else they may be too small to be seen
  ]
end



;;;;;;;;;; GO PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-presources
; Generates the amount of resources available in each patch (presources) depending on the quality
; of that patch (pquality). Patches resources are re-set to their value ± 10% at each tick,
; independently of resources level at the end of the previous tick. That is, resources do not
; build up (non-exploited resources decay) and are fully replenished each tick/season.
; Presources are determined by : - pquality (good, medium or bad patches) - the landscape mean
; amount of resources (i.e. slider-defined presources-growth) that allows setting whether the
; environment is rich or poor in resources. Note that the fact that presources are reset to
; their basis value ± 10 % at each tick introduces some spatial and temporal heterogeneity
; even in uniform landscapes

  ask patches
  [
    let base-resource-input (pquality * presources-growth)
      ; the basal resources of a patch are determined by pquality and the level of
      ; resource growth per tick
    set presources int ( base-resource-input + ((1 - random-float 2) * 0.1 *
        base-resource-input) )
          ; variance is +/- 0-10% of resource input
  ]
end


to partition-presources
; This function uses 3 reporters to calculate the competitive scores and size scores of all
; agents present in each patch and assign resources to them accordingly. All agents present
; are processed from smallest to largest. Each agent receives resources up to the limits
; prescribed by its relative competitive ability (dependent) and its size compared to a
; fully mature colony (independent) - whichever is lower. Both competitive ability and
; maturity based resource acquisition are defined by functions in reporters and the
; strength of these effects can be adjusted by changing values of constants in each.
; Larger agents are thus able to acquire resources that were not harvested by smaller
; agents because of size limited ability (independent) even if this exceeds their initial
; competitive share (but not their harvesting ability). Small agents cannot exploit all
; resources in a patch even if they are alone until they reach sufficient size.

  ask patches with [any? turtles-here]
    ; Identifies patches with turtles
  [
    let all-competitors turtles-here
      ; defines an agent set comprising all turtles present on the patch
    let remaining-competitors turtles-here
      ; defines a second agent set that will be modified along this procedure a
      ;turtles will be processed
    while [any? remaining-competitors]
      ; loops until all turtles have been processed (patch by patch)
    [
      let smallest-agent min-one-of remaining-competitors [tresources]
        ; finds the smallest turtle (random if tied for tresources)
      let share-of-presources int ( (competitiveness smallest-agent) /
        (summed-competitiveness all-competitors) * presources )
          ; determines the amount of resources it wins

      ask smallest-agent
      [
        ifelse (share-of-presources > tresources * max-growth-rate)
          ; if the amount of resources won exceeds the collection capacity of the agent
          ; (= tresources * max-growth-rate):
        [set tresources tresources + tresources * max-growth-rate]
          ; - the agent collects as much resources it can
        [set tresources tresources + share-of-presources]
          ; - otherwise it collects all it won
      ]
      set remaining-competitors remaining-competitors with [self != smallest-agent]
        ; removes processed agent (smallest agent) from agentset then loop
    ]
  ]
end



;;;;;;;;;; REPORTS OF partition-presources ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report competitiveness [x]
; this reports the competitive score of the focal agent, based on tresources relative to maturity
; threslhold. Large agents gain more when competitiveness is an exponential function of tresources
; than when it is a linear function, and they gain less when it is a logarithmic function.

  if presources-partitioning = "linear"      [report (         ([tresources] of x)
      / ([maturity-threshold] of x)  )]
  if presources-partitioning = "exponential" [report (exp (5 * ([tresources] of x)
        / ([maturity-threshold] of x)) )]
  if presources-partitioning = "logarithmic" [report (ln  (1 + ([tresources] of x)
        / ([maturity-threshold] of x)) )]
end


to-report summed-competitiveness [y]
; this reports the sum of competitive scores of all agents present in patch, including the
; smallest focal agent

  let cumulative-score 0
  ask y [set cumulative-score cumulative-score + competitiveness self]
  report cumulative-score
end



to maintain-turtle
; updates tresources by paying slider-defined maintenance cost

  ask turtles [set tresources floor (tresources * (1 - maintenance-cost / 100))]
end


to reproduce
; creates new turtles based on the reproductive strategy of the parent
; turtles allocate tresources (~body size) to reproduction once they reach reproductive size,
; according to the reproductive investment slider. Two strategies are implemented : minimum
; number of offspring and maximum number of offspring. Minimum creates a single offspring with
; body-size equal to all resources invested in reproduction. Maximum creates the most offspring
; possible of minimum body size (set by slider) and returns lefover resources to the parent
;
; NOTE: the dispersal procedure below includes three options: random dispersal, dispersal to
; patches with highest resources, dispersal to patch with lowest competition. The former two
; strategies requires agents to sample their environment (i.e. patches) before dispersing hence
; the procedure (1) finds a patch, then (2) the agent disperse to the centre of that patch
; and (3) it finishes dispersing within the patch. This means that the dispersing distance of
; any given agent is approximated and may be slightly over- or under-estimated. However, this
; has no effect on the results given that dispersal distances are drawn from a uniform
; distribution between 0 and "dispersal-range". When the dispersing agent cannot sample other
; patches it disperses randomly its dispersal-range in any direction.
;


  clear-links
    ; links are drawn between mother and offsprings to visualise dispersal.
    ; Here links of previous clicks are deleted for clarity

  ask turtles with [tresources > maturity-threshold]
    ; ask mature agents
  [
    let resources-for-reproduction tresources * (reproductive-investment / 100)
      ; calculate resources available for reproduction based on % investment
    let invested-resources 0
    let mother who
      ; records the unique ID number of the reproducing agent
    ifelse offspring-number = "min-nb-offspring"
      ; reproductive method for single offspring of large size
    [
      if random 100 > dispersal-risk
        ; only new agents surviving dispersal are plotted on the map. The others die during
        ; dispersal
      [
        hatch 1
        [
          set tresources resources-for-reproduction
            ; all available resources are invested in one single propagule
          set age 0
          ifelse patch-chosen dispersal-range = nobody
            ; if there is no patch to disperse to other than the current one
          [
            right random-float 360
              ; the agent takes a random direction
            forward random-float dispersal-range
              ; and disperse in that direction (i.e. within the patch)
          ]
          [
              ; if there is another patch to disperse to (chosen randomly or selectively)
            move-to patch-chosen dispersal-range
              ; the agent disperses to the center of that patch
            setxy (xcor - 0.5 + random-float 1) (ycor - 0.5 + random-float 1)
              ; and finishes dispersing within that patch
          ]
          if show-dispersal-s1? and breed = s1-agents [create-link-from turtle mother
            [set color ([color] of turtle mother)]]
          if show-dispersal-s2? and breed = s2-agents [create-link-from turtle mother
            [set color ([color] of turtle mother)]]
        ]
      ]
      set invested-resources resources-for-reproduction
        ; calculates resources the mother invested in her single large offspring
    ]
      ; end min-offspring ifelse

    [
        ; max offspring strategy produces the maximum number of offspring of minimum size
        ; with leftover resources returned to the parent
      let offspring-produced floor (resources-for-reproduction / min-offspring-resources)
      let a 1
      while [a <= offspring-produced]
        ; loop producing n offspring
      [
        if random 100 > dispersal-risk
          ; only new agents surviving dispersal are plotted on the map. The others die
        [
          hatch 1
          [
            set tresources min-offspring-resources
              ; each of the n offspring receives slider-defined resources
            set age 0
            ifelse patch-chosen dispersal-range = nobody
              ; if there is no patch to disperse to other than the current one
            [
              right random-float 360
                ; the agent takes a random direction
              forward random-float dispersal-range
                ; and disperses in that direction (i.e. within the patch)
            ]
            [
                ; if there is another patch to disperse to (chosen randomly or selectively)
              move-to patch-chosen dispersal-range
                ; the agent disperses to the center of that patch
              setxy (xcor - 0.5 + random-float 1) (ycor - 0.5 + random-float 1)
                ; and finishes dispersing within that patch
            ]
            if show-dispersal-s1? and breed = s1-agents [create-link-from turtle mother
              [set color ([color] of turtle mother)]]
            if show-dispersal-s2? and breed = s2-agents [create-link-from turtle mother
              [set color ([color] of turtle mother)]]
          ]
        ]
        set a a + 1
          ; loops until all new agents have been generated
        set invested-resources invested-resources + min-offspring-resources
          ; calculates resources the mother invested in her n offspring
          ; (including those that die during dispersal)
      ]
        ; end while
    ]
      ; end max offspring ifelse

    set tresources tresources - invested-resources
      ; calculate resources the mother has left
  ]

  ifelse show-turtle-size?
  [update-turtle-size]
  [ask turtles [set size 0.3]]
end


;;;;;;;;;; REPORTS OF reproduce ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report patch-chosen [potential-dispersal-range]
; reports the patch chosen to disperse to according to selected critera and range
; - "low-competition" = dispersal to patch with the smallest biomass (tresources) of agents
; - "high-resources" = dispersal to patch with most presources
; - "random" = random dispersal
; note: when range is lower than 1 reports NOBODY instead of reporting SELF

  if dispersal-patch-choice = "low-competition" [report min-one-of (patches in-radius potential-dispersal-range)
    [sum [tresources] of turtles-here]]
  if dispersal-patch-choice = "high-resources"  [report max-one-of (patches in-radius potential-dispersal-range)
    [presources]]
  if dispersal-patch-choice = "random"          [report     one-of (patches in-radius potential-dispersal-range)]
end



to live-or-die
  ask turtles
    ; turtles may die of old age or starvation
  [
    set age age + 1
    if age = death-age or tresources <= 0 [die]
  ]

  ask patches [if random 100 < disturbance [ask turtles-here [die]]]
    ; turtles may die from patch level stochastic events
end





to launch-invasion
; introduces an single 'invading' agent to an ongoing simulation at time t, where t is user
; defined (and should be when population is at equilibrium to test ESS conditions).
; the invading agent is automatically determined to be the strategy absent at the start of
; the simulation.
;

  if (start-agents-s1 = 0)
  [
    create-s1-agents 1
    [
      set tresources start-resources
      set age 1
      set death-age longevity
      set maintenance-cost maintenance
      set max-growth-rate max-growth
      set maturity-threshold maturity-age
      set reproductive-investment investment-in-reproduction
      set offspring-number offspring-number-s1
      set dispersal-range dispersal-range-s1
      set dispersal-risk dispersal-risk-s1
      set dispersal-patch-choice dispersal-patch-choice-s1
      setxy random-xcor random-ycor
      set color red
    ]
  ]

  if (start-agents-s2 = 0)
  [
    create-s2-agents 1
    [
      set tresources start-resources
      set age 1
      set death-age longevity
      set maintenance-cost maintenance
      set max-growth-rate max-growth
      set maturity-threshold maturity-age
      set reproductive-investment investment-in-reproduction
      set offspring-number offspring-number-s2
      set dispersal-range dispersal-range-s2
      set dispersal-risk dispersal-risk-s2
      set dispersal-patch-choice dispersal-patch-choice-s2
      setxy random-xcor random-ycor
      set color blue
    ]
  ]
end


to-report check-any-strategy-dead
  let x false

  ifelse not invasion?
  ; when not simulating an invasion, stops the simulation if a strategy that was initially
  ; present became eradicated
  [
    if (start-agents-s1 > 0 and count s1-agents = 0) or (start-agents-s2 > 0 and count
      s2-agents = 0) [set x true]
  ]

  ; when simulating an invasion, stops the simulation if a strategy is eradicated
  [
    if ticks > invasion-time
    [
      if count s1-agents = 0 or count s2-agents = 0 [set x true]
    ]
  ]

  report x
end
@#$#@#$#@
GRAPHICS-WINDOW
1431
10
2059
639
-1
-1
20.0
1
10
1
1
1
0
1
1
1
-15
15
-15
15
1
1
1
ticks
30.0

BUTTON
598
42
667
75
Go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
527
42
596
75
Setup
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

BUTTON
669
42
739
75
Go1
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
13
14
163
32
Strategy 1
14
0.0
1

TEXTBOX
262
10
412
28
Strategy 2
14
0.0
1

TEXTBOX
13
500
163
518
Environment
14
0.0
1

SLIDER
266
378
510
411
maturity-age
maturity-age
0
1000
700.0
1
1
tresources
HORIZONTAL

SLIDER
9
341
251
374
longevity
longevity
0
100
10.0
1
1
ticks
HORIZONTAL

SLIDER
9
376
252
409
maintenance
maintenance
0
100
5.0
1
1
% tresources
HORIZONTAL

SLIDER
266
342
511
375
investment-in-reproduction
investment-in-reproduction
0
100
50.0
1
1
%
HORIZONTAL

SLIDER
271
523
509
556
disturbance
disturbance
0
100
5.0
1
1
%
HORIZONTAL

SLIDER
267
451
511
484
min-offspring-resources
min-offspring-resources
0
500
10.0
1
1
tresources
HORIZONTAL

TEXTBOX
13
318
317
336
Strategy independent turtle variables
14
0.0
1

SLIDER
10
186
254
219
dispersal-range-s1
dispersal-range-s1
0
30
1.0
0.1
1
patches
HORIZONTAL

SLIDER
261
187
506
220
dispersal-range-s2
dispersal-range-s2
0
30
30.0
0.1
1
patches
HORIZONTAL

SLIDER
9
222
255
255
dispersal-risk-s1
dispersal-risk-s1
0
100
10.0
1
1
% mortality
HORIZONTAL

SLIDER
12
102
253
135
start-agents-s1
start-agents-s1
0
1000
100.0
1
1
NIL
HORIZONTAL

SLIDER
260
103
507
136
start-agents-s2
start-agents-s2
0
1000
100.0
1
1
NIL
HORIZONTAL

INPUTBOX
11
42
254
102
Label-s1
DCF
1
0
String

INPUTBOX
260
43
507
103
Label-s2
ICF
1
0
String

SLIDER
267
415
511
448
start-resources
start-resources
30
1000
300.0
1
1
tresources
HORIZONTAL

CHOOSER
11
139
253
184
offspring-number-s1
offspring-number-s1
"min-nb-offspring" "max-nb-offspring"
0

CHOOSER
261
140
507
185
offspring-number-s2
offspring-number-s2
"min-nb-offspring" "max-nb-offspring"
1

SWITCH
1582
665
1730
698
show-dispersal-s1?
show-dispersal-s1?
1
1
-1000

SWITCH
1583
700
1731
733
show-dispersal-s2?
show-dispersal-s2?
1
1
-1000

SLIDER
262
222
506
255
dispersal-risk-s2
dispersal-risk-s2
0
100
90.0
1
1
% mortality
HORIZONTAL

PLOT
761
10
1088
216
Number of agents
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"all" 1.0 0 -16777216 true "" "plot count turtles"
"s1" 1.0 0 -2674135 true "" "plot count s1-agents"
"s2" 1.0 0 -13345367 true "" "plot count s2-agents"

SLIDER
7
523
253
556
presources-growth
presources-growth
0
1000
300.0
1
1
NIL
HORIZONTAL

MONITOR
603
501
672
546
nb s1
count s1-agents
17
1
11

MONITOR
674
501
743
546
nb s2
count s2-agents
17
1
11

MONITOR
531
501
601
546
individuals
count turtles
17
1
11

PLOT
1096
10
1423
216
Mean resources (~body size) of agents
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"all" 1.0 0 -16777216 true "" "if any? turtles [plot mean [tresources] of turtles]"
"s1" 1.0 0 -2674135 true "" "ifelse not any? s1-agents [plot 0] [plot mean [tresources] of s1-agents]"
"s2" 1.0 0 -13345367 true "" "ifelse not any? s2-agents [plot 0] [plot mean [tresources] of s2-agents]"

SWITCH
1734
700
1900
733
show-presources-value?
show-presources-value?
1
1
-1000

PLOT
531
179
740
329
max patch resources
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot max [presources] of patches"

PLOT
530
337
739
487
mean patch resources
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [presources] of patches"

MONITOR
532
553
604
598
Ticks
ticks
17
1
11

CHOOSER
270
558
509
603
landscape-aggregation
landscape-aggregation
"uniform" "random" 0 0.5 1
1

SWITCH
1734
665
1902
698
show-pquality-value?
show-pquality-value?
1
1
-1000

TEXTBOX
1592
741
1917
759
pquality does not fluctuate over time, presources does
12
0.0
1

MONITOR
531
130
607
175
mean pquality
mean [pquality] of patches
2
1
11

SWITCH
530
86
657
119
stop-on-win?
stop-on-win?
0
1
-1000

SWITCH
1432
664
1579
697
show-turtle-size?
show-turtle-size?
0
1
-1000

PLOT
1095
219
1422
424
Mean age of agents
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"all" 1.0 0 -16777216 true "" "if any? turtles [plot mean [age] of turtles]"
"s1" 1.0 0 -2674135 true "" "ifelse not any? s1-agents [plot 0] [plot mean [age] of s1-agents]"
"s2" 1.0 0 -13345367 true "" "ifelse not any? s2-agents [plot 0] [plot mean [age] of s2-agents]"

TEXTBOX
1101
641
1422
693
Note that \"show-dispersal\", \"mean body size\" and \"mean age\" do not include the propagules that die during dispersal.
14
0.0
1

PLOT
762
428
1088
633
Number of agents per patch class
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"s1 pclass 1" 1.0 0 -2674135 true "" "ifelse not any? s1-agents [plot 0] [plot count s1-agents-on patches with [pclass = 1]]"
"s2 pclass 1" 1.0 0 -13345367 true "" "ifelse not any? s2-agents [plot 0] [plot count s2-agents-on patches with [pclass = 1]]"
"s1 pclass 2" 1.0 0 -1604481 true "" "ifelse not any? s1-agents [plot 0] [plot count s1-agents-on patches with [pclass = 2]]"
"s2 pclass 2" 1.0 0 -8020277 true "" "ifelse not any? s2-agents [plot 0] [plot count s2-agents-on patches with [pclass = 2]]"

PLOT
762
219
1087
424
Mean quality of occupied patches
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"s1" 1.0 0 -2674135 true "" "ifelse not any? s1-agents [plot 0] [plot mean [pquality] of patches with [count s1-agents-here > 0]]"
"s2" 1.0 0 -13345367 true "" "ifelse not any? s2-agents [plot 0] [plot mean [pquality] of patches with [count s2-agents-here > 0]]"

SLIDER
7
557
253
590
prop-class1-patches
prop-class1-patches
0
1
0.5
0.1
1
NIL
HORIZONTAL

SLIDER
7
592
251
625
class2-patches-relative-pquality
class2-patches-relative-pquality
0
1
0.5
0.05
1
NIL
HORIZONTAL

MONITOR
611
130
683
175
Good habitat %
count patches with [pquality = 1] / count patches
1
1
11

CHOOSER
9
257
256
302
dispersal-patch-choice-s1
dispersal-patch-choice-s1
"low-competition" "high-resources" "random"
2

CHOOSER
263
257
506
302
dispersal-patch-choice-s2
dispersal-patch-choice-s2
"low-competition" "high-resources" "random"
2

CHOOSER
8
448
253
493
presources-partitioning
presources-partitioning
"linear" "exponential" "logarithmic"
0

SLIDER
9
411
251
444
max-growth
max-growth
0
50
5.0
.5
1
x
HORIZONTAL

TEXTBOX
276
608
511
778
- 'uniform' yields a continuous landscape of medium patches.\n- 'random' and 'aggregated' yield landscapes of good and bad patches\n- Numbers refer to Hurst exponent values used to generate more or less  'aggregated' landscapes (0 = lowest aggregation, 1 = highest).
14
0.0
1

PLOT
1096
429
1423
633
Nb large agents per patch class
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"s1 pclass 1" 1.0 0 -2674135 true "" "ifelse not any? s1-agents [plot 0] [plot count (s1-agents-on patches with [pclass = 1]) with [tresources >= (maturity-age * 0.75)]]"
"s1 pclass 2" 1.0 0 -1069655 true "" "ifelse not any? s1-agents [plot 0] [plot count (s1-agents-on patches with [pclass = 2]) with [tresources >= (maturity-age * 0.75)]]"
"s2 pclass 1" 1.0 0 -14070903 true "" "ifelse not any? s2-agents [plot 0] [plot count (s2-agents-on patches with [pclass = 1]) with [tresources >= (maturity-age * 0.75)]]"
"s2 pclass 2" 1.0 0 -8020277 true "" "ifelse not any? s2-agents [plot 0] [plot count (s2-agents-on patches with [pclass = 2]) with [tresources >= (maturity-age * 0.75)]]"

SWITCH
8
638
108
671
invasion?
invasion?
1
1
-1000

INPUTBOX
112
638
194
698
invasion-time
500.0
1
0
Number

TEXTBOX
12
706
226
757
Invasion introduces one individual of the strategy absent at the start, at the tick entered
14
0.0
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

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

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

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
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

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
  <experiment name="invasion 1.5.0 base" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3000"/>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Label-s1">
      <value value="&quot;DCF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Label-s2">
      <value value="&quot;ICF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s1">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s2">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s1">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s2">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s2">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s1">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s1">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s2">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s1">
      <value value="&quot;min-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s2">
      <value value="&quot;max-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-offspring-resources">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-growth">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="longevity">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maintenance">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="investment-in-reproduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maturity-age">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-resources">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion-time">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-growth">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stochastic-mortality">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-aggregation">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-class1-patches">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="class1-patch-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="class2-patch-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="class2-patches-relative-pquality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-partitioning">
      <value value="&quot;linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-on-win?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s2?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-turtle-size?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s1?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-pquality-value?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-presources-value?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="invasion by DCF - longevity" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3000"/>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="invasion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion-time">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-growth">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stochastic-mortality">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="longevity">
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
      <value value="8"/>
      <value value="9"/>
      <value value="10"/>
      <value value="15"/>
      <value value="20"/>
      <value value="25"/>
      <value value="30"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Label-s1">
      <value value="&quot;DCF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Label-s2">
      <value value="&quot;ICF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s1">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s2">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s1">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s2">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s2">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s1">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s1">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s2">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s1">
      <value value="&quot;min-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s2">
      <value value="&quot;max-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-offspring-resources">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-growth">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="investment-in-reproduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maturity-age">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-resources">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maintenance">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-aggregation">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-class1-patches">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="class1-patch-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="class2-patch-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="class2-patches-relative-pquality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-partitioning">
      <value value="&quot;linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-on-win?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s2?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-turtle-size?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s1?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-pquality-value?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-presources-value?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="uniform ref run comp" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <metric>count s1-individuals</metric>
    <metric>count s2-individuals</metric>
    <metric>sum [tresources] of s1-individuals</metric>
    <metric>sum [tresources] of s2-individuals</metric>
    <metric>count s1-individuals-on patches with [pclass = 1]</metric>
    <metric>count s1-individuals-on patches with [pclass = 2]</metric>
    <metric>count s2-individuals-on patches with [pclass = 1]</metric>
    <metric>count s2-individuals-on patches with [pclass = 2]</metric>
    <metric>count (s1-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s1-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s2-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s2-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of s1-individuals-on patches with [pclass = 1]</metric>
    <metric>sum [tresources] of s1-individuals-on patches with [pclass = 2]</metric>
    <metric>sum [tresources] of s2-individuals-on patches with [pclass = 1]</metric>
    <metric>sum [tresources] of s2-individuals-on patches with [pclass = 2]</metric>
    <metric>sum [tresources] of (s1-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s1-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s2-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s2-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <enumeratedValueSet variable="maturity-age">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s1">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maintenance">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="investment-in-reproduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-aggregation">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s2?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion-time">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s2">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-pquality-value?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-partitioning">
      <value value="&quot;linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s1">
      <value value="&quot;min-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-turtle-size?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-growth">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-class1-patches">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-offspring-resources">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-resources">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-on-win?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s1">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-growth">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Label-s1">
      <value value="&quot;DCF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-presources-value?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s1?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s2">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s2">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="longevity">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s1">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s2">
      <value value="&quot;max-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Label-s2">
      <value value="&quot;ICF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s1">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="class2-patches-relative-pquality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s2">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="100% bad" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <metric>count s1-individuals</metric>
    <metric>count s2-individuals</metric>
    <metric>sum [tresources] of s1-individuals</metric>
    <metric>sum [tresources] of s2-individuals</metric>
    <metric>count s1-individuals-on patches with [pclass = 1]</metric>
    <metric>count s1-individuals-on patches with [pclass = 2]</metric>
    <metric>count s2-individuals-on patches with [pclass = 1]</metric>
    <metric>count s2-individuals-on patches with [pclass = 2]</metric>
    <metric>count (s1-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s1-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s2-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s2-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of s1-individuals-on patches with [pclass = 1]</metric>
    <metric>sum [tresources] of s1-individuals-on patches with [pclass = 2]</metric>
    <metric>sum [tresources] of s2-individuals-on patches with [pclass = 1]</metric>
    <metric>sum [tresources] of s2-individuals-on patches with [pclass = 2]</metric>
    <metric>sum [tresources] of (s1-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s1-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s2-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s2-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <enumeratedValueSet variable="presources-growth">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-class1-patches">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s1">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-growth">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s1">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-aggregation">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion-time">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s1">
      <value value="&quot;min-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s1?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s2">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="class2-patches-relative-pquality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="longevity">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-presources-value?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maintenance">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Label-s1">
      <value value="&quot;DCF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-on-win?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-offspring-resources">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="investment-in-reproduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s2">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-resources">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s2?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Label-s2">
      <value value="&quot;ICF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s2">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s1">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-partitioning">
      <value value="&quot;linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s2">
      <value value="&quot;max-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s1">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-pquality-value?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maturity-age">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-turtle-size?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s2">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="resources uniform rerun - ref" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <metric>count s1-individuals</metric>
    <metric>count s2-individuals</metric>
    <metric>sum [tresources] of s1-individuals</metric>
    <metric>sum [tresources] of s2-individuals</metric>
    <metric>count s1-individuals-on patches with [pclass = 1]</metric>
    <metric>count s1-individuals-on patches with [pclass = 2]</metric>
    <metric>count s2-individuals-on patches with [pclass = 1]</metric>
    <metric>count s2-individuals-on patches with [pclass = 2]</metric>
    <metric>count (s1-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s1-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s2-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s2-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of s1-individuals-on patches with [pclass = 1]</metric>
    <metric>sum [tresources] of s1-individuals-on patches with [pclass = 2]</metric>
    <metric>sum [tresources] of s2-individuals-on patches with [pclass = 1]</metric>
    <metric>sum [tresources] of s2-individuals-on patches with [pclass = 2]</metric>
    <metric>sum [tresources] of (s1-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s1-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s2-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s2-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <enumeratedValueSet variable="dispersal-risk-s1">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-presources-value?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-pquality-value?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s1?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s2">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-turtle-size?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-on-win?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="class2-patches-relative-pquality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="investment-in-reproduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Label-s2">
      <value value="&quot;ICF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-growth">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-class1-patches">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Label-s1">
      <value value="&quot;DCF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s2">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-aggregation">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s2">
      <value value="&quot;max-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-offspring-resources">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-resources">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maturity-age">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s2?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maintenance">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion-time">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="longevity">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s1">
      <value value="&quot;min-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s1">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-growth">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s1">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s2">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s2">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s1">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-partitioning">
      <value value="&quot;linear&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="invasion 1.5.1 ICF disturbance" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3000"/>
    <metric>count s1-individuals</metric>
    <metric>count s2-individuals</metric>
    <metric>sum [tresources] of s1-individuals</metric>
    <metric>sum [tresources] of s2-individuals</metric>
    <metric>count s1-individuals-on patches with [pclass = 1]</metric>
    <metric>count s1-individuals-on patches with [pclass = 2]</metric>
    <metric>count s2-individuals-on patches with [pclass = 1]</metric>
    <metric>count s2-individuals-on patches with [pclass = 2]</metric>
    <metric>count (s1-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s1-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s2-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>count (s2-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of s1-individuals-on patches with [pclass = 1]</metric>
    <metric>sum [tresources] of s1-individuals-on patches with [pclass = 2]</metric>
    <metric>sum [tresources] of s2-individuals-on patches with [pclass = 1]</metric>
    <metric>sum [tresources] of s2-individuals-on patches with [pclass = 2]</metric>
    <metric>sum [tresources] of (s1-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s1-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s2-individuals-on patches with [pclass = 1]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <metric>sum [tresources] of (s2-individuals-on patches with [pclass = 2]) with [tresources &gt;= (maturity-age * 0.75)]</metric>
    <enumeratedValueSet variable="Label-s1">
      <value value="&quot;DCF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Label-s2">
      <value value="&quot;ICF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s1">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-individuals-s2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s1">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-range-s2">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s2">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-patch-choice-s1">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s1">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal-risk-s2">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s1">
      <value value="&quot;min-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offspring-number-s2">
      <value value="&quot;max-nb-offspring&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-offspring-resources">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-growth">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="longevity">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maintenance">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="investment-in-reproduction">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maturity-age">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-resources">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion-time">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-growth">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disturbance">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="7"/>
      <value value="10"/>
      <value value="15"/>
      <value value="20"/>
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-aggregation">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-class1-patches">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="class2-patches-relative-pquality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presources-partitioning">
      <value value="&quot;linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-on-win?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s2?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-turtle-size?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-dispersal-s1?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-pquality-value?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-presources-value?">
      <value value="false"/>
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
