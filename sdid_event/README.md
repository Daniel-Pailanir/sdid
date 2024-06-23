# sdid_event

A Stata module to implement event study analysis with `sdid`.

# Setup

```stata
net install sdid_event, from("https://raw.githubusercontent.com/DiegoCiccia/sdid/main/sdid_event") replace

```

# Example

```stata
clear
local GG = 19
local TT = 20
set seed 0
set obs `=`GG' * `TT''

qui do "sdid_event.ado"
gen G = mod(_n-1,`GG') + 1
bys G: gen T = _n

gen D = T > mod(G, 4) + 1 & G >= `GG'/4
gen Y = uniform() * (1 + D + 10*D*T)
```

Basic syntax - all the effects are retrieved:

```stata
sdid_event Y G T D
```

**effects()** - limiting the number of reported estimates:

```stata
sdid_event Y G T D, effects(5)
```

**disag** - reporting cohort-specific dynamic treatment effect estimates:

```stata
sdid_event Y G T D, effects(5) disag
```


