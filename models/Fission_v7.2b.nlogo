;; initial_colony located at the centre, with all resources except foragers
;; potential nesting sites = nests
;; new_colonies formed by fission of initial_colony and located at chosen_nests
;;
;;
;; outputs:
;;
;; 7 variables = number_new_colonies; mean_resources_new_colonies (mean_resources_chosen_nests) ; CV (coefficient of variation between new_colonies)
;; mean_quality_chosen_nests ; mean_quality_all_nests ; mean_distance_chosen_nests ; mean_distance_all_nests
;;
;; 3 lists = raw_resources_new_colonies; raw_quality_chosen_nests; raw_distance_chosen_nests
;; the 3 lists are ranked by decreasing resources allocated to new_colonies (= nests receiving resources, and nests that received no resources are not listed)
;;
;;


;;;;;;;;;; DEFINITION OF VARIABLES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

globals
[
  number_foragers

  ;; for raw outputs
  number_new_colonies
  chosen_nests                          ; temp agentset to rank chosen_nests by size
  raw_resources_new_colonies
  raw_quality_all_nests
  raw_quality_chosen_nests
  raw_distance_chosen_nests

  ;; for "aggregated" outputs
  mean_resources_new_colonies
  CV                                     ; coefficient of variation in resources between chosen_nests
  mean_quality_chosen_nests
  mean_quality_all_nests
  mean_distance_chosen_nests
  mean_distance_all_nests

  visiting_workers_all_nests
  visited_nests_all_workers
]

breed [initial_colonies initial_colony]  ; initial colony, from where foragers will transport resources to new colonies founded by colony fission
breed [nests nest]                       ; these are the potential nesting sites, where new_colonies could be founded
breed [workers worker]                   ; these are the fraction of total workers that are actually responsible for colony fission is resource allocation to new nest(s)
breed [end_workers end_worker]           ; these are workers once they self allocated to a new nest

turtles-own [resources]                  ; initial_colony, nests and foragers have or transport resources (ie other workers, queen and brood)
nests-own
[
  quality                      ; nests may vary in quality
  visiting_workers
]
workers-own
[
  state                                  ; workers are in 1 of 4 states: exploring, returning to the initial colony, transporting resource to a chosen nest, self_allocating oneself to a chosen nest
  motivated_to_transport?
  current_visited_nest                   ; nest visited last
  visited_nests                          ; agentset of visited nests
  chosen_nest                            ; ID of chosen nest
  assessed_quality_of_chosen_nest
]

end_workers-own [visited_nests]

;;;;;;;;;; SETUP ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all   ; to remove/annotate when running wth OpenMole
  clear-output
  setup-patches
  setup-initial_colony
  setup-nests
  setup-workers
  reset-ticks
end



;;;;;;;;;; GO ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go
  move                                                 ; workers move to (1) explore the map, (2) return to the initial colony and (3) transport resources to a chosen nest

  if show_nests_size? [update_nests_size]              ; option to graphically show resource build up in chosen nests (and depletion in initial colony)

  tick

  if not any? workers or ticks = max_ticks
  [
    output_data_to_interface
    if write_data? [write_data_to_file]
    stop
  ]
end


;;;;;;;;;; SETUP PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-patches
  ask patches [set pcolor green]
end


to setup-initial_colony                                  ; sets up the initial_colony in the centre of the map with all resources
  set number_foragers int (size_initial_colony * percentage_foragers / 100)

  create-initial_colonies 1
  [
    setxy 0 0
    set shape "circle"
    set color black
    set resources size_initial_colony - number_foragers
    set label resources
    if show_nests_size? [update_nests_size]
  ]
end


to setup-nests
  ifelse position_nests = "equidistant"
  [                                                                              ; sets up nests in a circle around the initial_colony, at a distance of 25 patches
    create-ordered-nests number_nests [jump 25]
  ]
  ; end of ifelse position = equidistant

  [                                                                              ; sets up nests in random position BUT not in an exclusion radius around the intial nest
    ask [patches in-radius exclusion_radius] of patch 0 0 [set pcolor lime]      ; visualises the exclusion radius (approximate as it show patches not actual exclusion radius)
    create-nests number_nests
    [
      setxy random-xcor random-ycor
      while [any? nests with [distance initial_colony 0 < exclusion_radius]]     ; nests generated in the exclusion radius change position until no nests are in the exclusion radius
      [
        ask nests with [distance initial_colony 0 < exclusion_radius]
        [
          setxy random-xcor random-ycor
          ; note that some nests may be located very close toone another. If one wishes, one could avoid this by re-creating nests falling within radius eg 1 of an existing nest
        ]
      ]
    ]
  ]
  ; end of ifelse position = random


  ask nests
  [
    set shape "circle"
    set color brown
    set quality int random-normal nests_quality nests_quality_SD
    if quality < 0 [set quality 0]
    set label quality
    set resources 0
    set visiting_workers no-turtles
    if show_nests_size? [update_nests_size]
  ]
end


to setup-workers
  create-workers number_foragers                                ; foragers are black bugs and start in the initial_colony
  [
    set shape "bug"
    setxy 0 0
    set color black
    set current_visited_nest nobody                             ; empty agent
    set visited_nests no-turtles                                ; empty agentset
    set chosen_nest nobody
    set state "exploring"
    set motivated_to_transport? false                           ; at the beginning, foragers are not motivated to transport but only to explore
  ]
end



;;;;;;;;;; GO PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to move
  ask workers
  [
    ifelse ticks < exploring_phase
    [
      explore
      chose_nests
    ]
    ; end of ifelse < exploring_phase

    [
      if ticks = exploring_phase and chosen_nest != nobody [set state "returning"]
      if state = "exploring" [explore chose_nests]
      if state = "returning" [return]
      if state = "transporting" [transport]
      if state = "self_allocating" [self_allocate]
    ]
    ; end of ifelse > exploring_phase
  ]
end

to explore                          ; workers explore the map with random walk with + or - 45° deviation
  ifelse random 100 < 50
  [right random 45]
  [left random 45]
  forward 0.75                      ; the lenght of movement must be less than the detection radius around nests, else an agent could "jump" over a nest without detecting it

  ; if the intial colony is depleted of resources the foragers who failed to find and chose a nest self-allocate themselves to one of the nests chosen by other foragers
  ; It is implicit that they received this information from informed foragers
  ; It could be modified so that they would go to best nests if foragers who found better nest are more inclined to inform others
  if [resources] of initial_colony 0 = 0 and chosen_nest = nobody
  [
    set chosen_nest one-of nests with [resources > 0]
    set state "self_allocating"
    set color yellow
  ]
end


to chose_nests
  if any? nests in-radius 1                                                      ; any worker within a radius 1 around a nest detects this nest
  [
    set current_visited_nest min-one-of nests in-radius 1 [distance self]        ; if a worker detects two nests she chooses the closest one

    ask current_visited_nest
    [
      set color blue
      set visiting_workers (turtle-set visiting_workers myself)
    ]

    set color blue
    set visited_nests (turtle-set visited_nests current_visited_nest)

    if choice_strategy = "best_nest" [set chosen_nest max-one-of visited_nests [int random-normal (quality) (quality * nest_quality_assessment_error / 100)] ]
    ; Note that the quality of nests that have been visited multiple times is not better assessed!
    if choice_strategy = "last_nest" [set chosen_nest current_visited_nest]
    if choice_strategy = "rotating_choice" [set chosen_nest one-of visited_nests]

    set current_visited_nest nobody

    if ticks > exploring_phase and chosen_nest != nobody [set state "returning"]
  ]
end


to return
  face initial_colony 0
  forward 0.75

  if distance initial_colony 0 < 1              ; if the initial_colony is at less than 1 unit of distance the worker is at the initial_colony
  [
    ifelse [resources] of initial_colony 0 != 0

    ; if the initial_colony is not empty of resources, the forager collects resources and starts going back to the chosen nest
    [
      ; updates rotating_choice at each return to the initial colony
      if choice_strategy = "rotating_choice" [set chosen_nest one-of visited_nests]

      check_motivation_to_transport
      if motivated_to_transport?
      [
        collect_resource
        set state "transporting"
      ]
    ]

    ; if the initial_colony is empty, the forager starts going back to the chosen nest in order to join it
    [
      set state "self_allocating"
      set color yellow
    ]
  ]
end



to check_motivation_to_transport
  if probability_of_transporting = "unconditionnal" [set motivated_to_transport? true]

  ; foragers that have chosen a nest transport ressources to it with a probability that depend on its quality
  ; => it's probabilistic if the quality of the chosen nest < 100, it's certain if quality >= 100
  if probability_of_transporting = "linear"
  [
    if assessed_quality_of_chosen_nest >= random 100 [set motivated_to_transport? true]
  ]

  if probability_of_transporting = "exponential"
  [
    if exp (assessed_quality_of_chosen_nest * exponential_factor) >= random 100 [set motivated_to_transport? true]
  ]

  if probability_of_transporting = "logarithmic"
  [
    if ln (1 + assessed_quality_of_chosen_nest) * logarithmic_factor >= random 100  [set motivated_to_transport? true]
  ]

  if probability_of_transporting = "logistic"
  [
    if 100 / (1 + exp (- (assessed_quality_of_chosen_nest - logistic_factor_mu) / logistic_factor_s)) >= random 100  [set motivated_to_transport? true]
  ]
end


to transport
  face chosen_nest
  forward 0.75
  if distance chosen_nest < 1                 ; if the chosen nest is at less than 1 unit of distance the worker is at the chosen nest
  [
    deposit_resource
    set motivated_to_transport? false         ; the motivation to transport is reset to false so that it is tested at each return to the initial colony
  ]
end


to self_allocate
  face chosen_nest
  forward 0.75
  if distance chosen_nest < 1                                    ; if the chosen nest is at less than 1 unit of distance the worker joins this nest
  [
    ask chosen_nest
    [
      set color yellow
      set resources (resources + 1)
      set label resources
    ]
    set breed end_workers
  ]
end


to collect_resource
  if [resources] of initial_colony 0 > 0
  [
    set resources 1
    set state "transporting"
    set color red

    ask initial_colony 0
    [
      set resources (resources - 1)
      set label resources
      if resources = 0
      [
        ask initial_colony 0 [set color yellow]
      ]
    ]
  ]
end


to deposit_resource
  set resources 0
  set state "returning"
  set color blue

  ask chosen_nest
  [
    set color red
    set resources (resources + 1)
    set label resources
  ]
end


to update_nests_size                                         ; scales size of nests proportional to their resources with max shown size of 3
  ask initial_colonies
  [
    set size 3 * resources / size_initial_colony
    if size < 1 [set size 1]                                 ; very small nests are shown larger than proportionally else they may be too small to be seen
  ]

  ask nests
  [
    set size 20 * resources / size_initial_colony
    if size < 1 [set size 1]
  ]
end


to output_data_to_interface
  ; this is useful for OpenMole to easily gets the outputs

  ; calculates outputs
  set number_new_colonies count nests with [resources > 0]
  set mean_resources_new_colonies mean [resources] of nests with [resources > 0]

  ifelse count nests with [resources > 0] >= 2
  [set CV (standard-deviation [resources] of nests with [resources > 0]) / (mean [resources] of nests with [resources > 0]) * 100]
  [set CV 0]

  set mean_quality_chosen_nests mean [quality] of nests with [resources > 0]
  set mean_quality_all_nests mean [quality] of nests
  set mean_distance_chosen_nests mean [distancexy 0 0] of nests with [resources > 0]
  set mean_distance_all_nests mean [distancexy 0 0] of nests

  set chosen_nests nests with [resources > 0]
  set raw_resources_new_colonies []
  set raw_quality_all_nests []
  set raw_quality_chosen_nests []
  set raw_distance_chosen_nests []

  set visiting_workers_all_nests []
  set visited_nests_all_workers []

  let i 0
  while [i <= number_nests]
  [
    ask nests with [who = i]
    [
      set raw_quality_all_nests lput quality raw_quality_all_nests
      ifelse any? visiting_workers
      [set visiting_workers_all_nests lput count visiting_workers visiting_workers_all_nests]
      [set visiting_workers_all_nests lput 0 visiting_workers_all_nests]
    ]
    set i i + 1
  ]

  let j 0
  while [j <= count end_workers]
  [
    ask end_workers with [who = j + number_nests]
    [
      ifelse any? visited_nests
      [set visited_nests_all_workers lput count visited_nests visited_nests_all_workers]
      [set visited_nests_all_workers lput 0 visited_nests_all_workers]
    ]
    set j j + 1
  ]

  while [any? chosen_nests]
  [
    ask max-one-of chosen_nests [resources]
    [
      set raw_resources_new_colonies lput resources raw_resources_new_colonies
      set raw_quality_chosen_nests lput quality raw_quality_chosen_nests
      set raw_distance_chosen_nests lput precision (distancexy 0 0) 3 raw_distance_chosen_nests
      set chosen_nests other chosen_nests
    ]
  ]

  ; writes outputs on interface
  output-type "number_new_colonies          " output-print number_new_colonies
  output-type "mean_resources_new_colonies  " output-print precision mean_resources_new_colonies 3
  output-type "CV                           " output-print precision CV 3
  output-type "mean_quality_chosen_nests    " output-print precision mean_quality_chosen_nests 3
  output-type "mean_quality_all_nests       " output-print precision mean_quality_all_nests 3
  output-type "mean_distance_chosen_nests   " output-print precision mean_distance_chosen_nests 3
  output-type "mean_distance_all_nests      " output-print precision mean_distance_all_nests 3
  output-type "raw_resources_new_colonies   " output-print raw_resources_new_colonies
  output-type "raw_quality_chosen_nests     " output-print raw_quality_chosen_nests
  output-type "raw_quality_all_nests        " output-print raw_quality_all_nests
  output-type "raw_distance_chosen_nests    " output-print raw_distance_chosen_nests
  output-type "visiting_workers_all_nests   " output-print visiting_workers_all_nests
  output-type "visited_nests_all_workers    " output-print visited_nests_all_workers

  ; use export-output ?
end


to write_data_to_file
  ; writes the data in two files:
  ; file_name_raw_data.csv lists raw measures (who, resources, quality, distance) for all chosen nests in rows
  ; file_name_mean_data.csv lists "synthetic" measures (number of chosen_nests, their mean resources, qualities and distances, and the coefficient of variation of their resources) on rows


  let temp_raw word "Raw_data_" file_name
  let raw_data_file_name word temp_raw ".csv"

  let temp_mean word "Mean_data_" file_name
  let mean_data_file_name word temp_mean ".csv"

  file-open raw_data_file_name

  let i 0
  file-type file_name file-type ";" file-type word "run_number_" behaviorspace-run-number file-type ";" file-type "nest_ID" file-type ";"
  while [i <= number_nests]
  [
    ask nests with [who = i] [file-type who file-type ";"]
    set i i + 1
  ]
  file-print ";"

  let j 0
  file-type file_name file-type ";" file-type word "run_number_" behaviorspace-run-number file-type ";" file-type "size_of_nest" file-type ";"
  while [j <= number_nests]
  [
    ask nests with [who = j] [file-type resources file-type ";"]
    set j j + 1
  ]
  file-print ";"

  let k 0
  file-type file_name file-type ";" file-type word "run_number_" behaviorspace-run-number file-type ";" file-type "quality_of_nest" file-type ";"
  while [k <= number_nests]
  [
    ask nests with [who = k] [file-type quality file-type ";"]
    set k k + 1
  ]
  file-print ";"

  let l 0
  file-type file_name file-type ";" file-type word "run_number_" behaviorspace-run-number file-type ";" file-type "distance_nest_to_initial_colony" file-type ";"
  while [l <= number_nests]
  [
    ask nests with [who = l] [file-type distancexy 0 0 file-type ";"]
    set l l + 1
  ]
  file-print ";"

  file-close


  file-open mean_data_file_name

  ; writes data on 2 lines : title above, and actual data below
  ; writes title of columns
  file-type "treatment" file-type ";"
  file-type "run_number" file-type ";"
  file-type "number_new_colonies" file-type ";"
  file-type "mean_resources_new_colonies" file-type ";"
  file-type "coefficient_of_variation" file-type ";"
  file-type "mean_quality_chosen_nests" file-type ";"
  file-type "mean_quality_all_nests" file-type ";"
  file-type "mean_distance_chosen_nests" file-type ";"
  file-type "mean_distance_all_nests" file-type ";"
  file-type "size_initial_colony" file-type ";"
  file-type "number_foragers" file-type ";"
  file-type "number_nests" file-type ";"
  file-type "nests_quality" file-type ";"
  file-type "nests_quality_SD" file-type ";"
  file-type "position_nests" file-type ";"
  file-type "exclusion_radius" file-type ";"
  file-type "nest_quality_assessment_error" file-type ";"
  file-type "choice_strategy" file-type ";"
  file-print "exploring_phase"

  ; writes data generated and values used for variables
  file-type file_name file-type ";"
  file-type word "run_number_" behaviorspace-run-number file-type ";"
  file-type number_new_colonies file-type ";"
  file-type mean_resources_new_colonies file-type ";"
  file-type CV file-type ";"
  file-type mean_quality_chosen_nests file-type ";"
  file-type mean_quality_all_nests file-type ";"
  file-type mean_distance_chosen_nests file-type ";"
  file-type mean_distance_all_nests file-type ";"
  file-type size_initial_colony file-type ";"
  file-type number_foragers file-type ";"
  file-type number_nests file-type ";"
  file-type nests_quality file-type ";"
  file-type nests_quality_SD file-type ";"
  file-type position_nests file-type ";"
  file-type exclusion_radius file-type ";"
  file-type nest_quality_assessment_error file-type ";"
  file-type choice_strategy file-type ";"
  file-print exploring_phase

  file-close
end
@#$#@#$#@
GRAPHICS-WINDOW
975
10
1593
629
-1
-1
10.0
1
10
1
1
1
0
0
0
1
-30
30
-30
30
1
1
1
ticks
30.0

BUTTON
257
45
322
78
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
257
10
323
43
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
257
81
323
114
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

MONITOR
661
10
731
55
Ticks
ticks
17
1
11

SLIDER
-1
10
242
43
size_initial_colony
size_initial_colony
0
10000
1000.0
1
1
NIL
HORIZONTAL

SLIDER
-1
46
243
79
percentage_foragers
percentage_foragers
0
100
10.0
1
1
%
HORIZONTAL

SLIDER
0
432
248
465
max_ticks
max_ticks
0
1000000
100000.0
1
1
ticks
HORIZONTAL

SLIDER
0
97
245
130
number_nests
number_nests
0
100
40.0
1
1
NIL
HORIZONTAL

CHOOSER
0
202
244
247
position_nests
position_nests
"random" "equidistant"
0

SLIDER
0
132
243
165
nests_quality
nests_quality
0
100
70.0
1
1
NIL
HORIZONTAL

SLIDER
0
248
246
281
exclusion_radius
exclusion_radius
0
30
10.0
1
1
patches
HORIZONTAL

SLIDER
0
167
243
200
nests_quality_SD
nests_quality_SD
0
100
50.0
1
1
NIL
HORIZONTAL

SWITCH
0
492
247
525
show_nests_size?
show_nests_size?
0
1
-1000

PLOT
257
268
611
418
Resources: temporal evolution
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
"Initial_colony" 1.0 0 -16777216 true "" "plot [resources] of initial_colony 0"
"Foragers" 1.0 0 -13345367 true "" "plot count workers"
"Sum in new_colonies" 1.0 0 -2674135 true "" "plot sum [resources] of nests"
"Largest new_colony" 1.0 0 -955883 true "" "plot max [resources] of nests"

PLOT
257
117
611
267
Workers status
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
"Exploring" 1.0 0 -16777216 true "" "plot count workers with [state = \"exploring\"]"
"Returning" 1.0 0 -13345367 true "" "plot count workers with [state = \"returning\"]"
"Transporting" 1.0 0 -2674135 true "" "plot count workers with [state = \"transporting\"]"
"Self_allocating" 1.0 0 -1184463 true "" "plot count workers with [state = \"self_allocating\"]"

SLIDER
0
397
249
430
exploring_phase
exploring_phase
0
10000
2000.0
1
1
ticks
HORIZONTAL

SLIDER
0
301
246
334
nest_quality_assessment_error
nest_quality_assessment_error
0
100
3.0
1
1
%
HORIZONTAL

CHOOSER
0
337
247
382
choice_strategy
choice_strategy
"best_nest" "last_nest" "rotating_choice"
0

PLOT
618
117
969
267
Resources = f (quality)
quality
resources
0.0
10.0
0.0
10.0
true
true
"set-plot-x-range 0 ceiling (1.1 * max [quality] of nests)" "clear-plot\nset-plot-x-range 0 ceiling (1.1 * max [quality] of nests)"
PENS
"Chosen nests" 1.0 2 -2674135 true "" "ask nests with [resources > 0] [plotxy quality resources]"
"Ignored nests" 1.0 2 -16777216 true "" "ask nests with [resources = 0] [plotxy quality resources]"

PLOT
617
268
968
418
Resources =f (distance)
distance
resources
0.0
10.0
0.0
10.0
true
true
"set-plot-x-range 0 ceiling (1.1 * max [distancexy 0 0] of nests)" "clear-plot\nset-plot-x-range 0 ceiling (1.1 * max [distancexy 0 0] of nests)"
PENS
"Chosen nests" 1.0 2 -2674135 true "" "ask nests with [resources > 0] [plotxy distancexy 0 0 resources]"
"Ignored nests" 1.0 2 -16777216 true "" "ask nests with [resources = 0] [plotxy distancexy 0 0 resources]"

MONITOR
816
10
966
55
Resources at initial_colony
[resources] of initial_colony 0
17
1
11

MONITOR
818
62
964
107
Resources in nests
sum [resources] of nests
17
1
11

MONITOR
662
63
809
108
Resources on workers
sum [resources] of workers
17
1
11

INPUTBOX
327
54
599
114
file_name
workers_and_brood
1
0
String

TEXTBOX
471
10
626
47
Write the name of the file where to save the data
14
15.0
1

SWITCH
328
10
462
43
write_data?
write_data?
1
1
-1000

OUTPUT
456
421
969
642
12

MONITOR
739
10
808
55
workers
count workers
17
1
11

CHOOSER
259
420
450
465
probability_of_transporting
probability_of_transporting
"unconditionnal" "linear" "exponential" "logarithmic" "logistic"
0

SLIDER
259
467
451
500
exponential_factor
exponential_factor
0
.3
0.05
0.005
1
NIL
HORIZONTAL

SLIDER
259
501
451
534
logarithmic_factor
logarithmic_factor
0
100
20.0
1
1
NIL
HORIZONTAL

SLIDER
259
538
451
571
logistic_factor_mu
logistic_factor_mu
0
100
80.0
1
1
NIL
HORIZONTAL

SLIDER
260
575
450
608
logistic_factor_s
logistic_factor_s
0
100
5.0
1
1
NIL
HORIZONTAL

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
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Resource allocation" repetitions="500" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <final>write_data_to_file</final>
    <exitCondition>sum [resources] of new_nests = workers_and_brood or ticks = max_ticks</exitCondition>
    <enumeratedValueSet variable="workers_and_brood">
      <value value="800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragers">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb_new_nests">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new_nests_quality">
      <value value="70"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new_nests_quality_SD">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Position_new_nests">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exclusion_radius">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nest_quality_assessment_error">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choice_strategy">
      <value value="&quot;best_new_nest&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exploring_phase">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_ticks">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show_nests_size?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Write_data?">
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
