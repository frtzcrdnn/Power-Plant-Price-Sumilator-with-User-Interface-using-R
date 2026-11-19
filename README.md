# Luzon Grid Dispatch and Price Simulator
 
An R Shiny app that finds the cheapest way to run the Luzon power plants for one day, then reads the wholesale electricity price of each 4-hour block out of the same solution. The optimization is a linear program solved with a Simplex implementation written from scratch in base R, with no LP library involved.
 
A full user manual in English and Filipino is in [`docs/150-UserManual-Cerdeña.pdf`](docs/150-UserManual-Cerdeña.pdf).
 
## What it does
 
The day is split into six 4-hour blocks, B1 (12-4 AM) through B6 (8 PM-12 MN). You pick which plant groups are available and, optionally, a daily CO₂ limit. The app then returns:
 
- the minimum total cost of the day, in pesos
- the output in MW of every plant group in every block
- the wholesale price of electricity in each block, in ₱/MWh
- every Simplex iteration, which you can step through in a second tab
## Requirements
 
- R 4.4 or newer. The app uses the `%||%` operator, which base R added in 4.4.0.
- The `shiny` package.
- RStudio is optional but makes launching a single keystroke.
```r
install.packages("shiny")
```
 
## Folder layout
 
The app reads its data and styling from relative paths, so the folders must look like this and R's working directory must be the project root:
 
```
.
├── CerdeñaProject.R
├── data
│   └── Primal.csv
├── ui
│   ├── custom_table.html
│   ├── header.html
│   └── styles.css
└── docs
    └── 150-UserManual-Cerdeña.pdf
```
 
## Running the app
 
In RStudio, open `CerdeñaProject.R` and press **Ctrl + Shift + S** (Source). The app opens in your browser.
 
From any R console instead:
 
```r
setwd("path/to/this/project")
source("CerdeñaProject.R")
shiny::runApp(app)
```
 
## Using it
 
1. **Plant groups.** Tick the groups that are available, or press **Check All**. **Reset** clears the selection, the CO₂ cap and any results.
2. **CO₂ cap (optional).** Tick *Limit daily CO₂ emissions* and enter a limit in tonnes per day. The default is 165,000.
3. **Solve.** Press **SOLVE**. The app switches to the Results tab.
The three **Test** buttons fill in the selection and CO₂ settings for the cases below. They do not solve anything by themselves, so press **SOLVE** after loading one.
 
### Results tab
 
| Block | Contents |
| --- | --- |
| Summary cards | Total minimum cost Z, number of Simplex iterations, number of unknowns, number of constraints |
| Output of each plant group | MW per group for B1 to B6, with row shading proportional to output and a column total |
| Price of electricity | Shadow price S<sub>t</sub> and the resulting price S<sub>t</sub> / 4 in ₱/MWh for each block |
| Shadow prices and check | Every non-zero shadow price with its constraint, RHS value b<sub>i</sub> and the product b<sub>i</sub> · S<sub>i</sub> |
 
The sum of the b<sub>i</sub> · S<sub>i</sub> column equals Z, which is the strong duality check.
 
### Simplex Tableau tab
 
Use **First**, **Prev**, **Next**, **Last**, the slider, or the play button under the slider to move through the iterations. A banner above the table names the pivot column, pivot row and pivot element for the current step. In the table, shadow price columns have brown headers, slack and unknown columns have navy headers, the pivot row is yellow, the pivot column is light blue and the pivot element is gold. A second table shows the basic solution, which is the bottom row of the current tableau.
 
## The model
 
Everything the solver needs lives in `data/Primal.csv`: 66 constraint rows, 42 decision variables and one objective row. The variables are the output x<sub>gt</sub> in MW of each plant group g in each block t, named like `Coal_B1` or `NatGas_B4`.
 
**Objective.** Minimize total generation cost. Each coefficient is the cost per MWh multiplied by the 4 hours in a block.
 
| Group | Nameplate capacity | Cost per MWh in the data | CO₂ coefficient (t per MW per block) |
| --- | --- | --- | --- |
| Coal | 7,000 MW | ₱3,800 to ₱3,990 | 3.8 |
| Natural Gas | 3,300 MW | ₱6,500 to ₱6,825 | 1.6 |
| Geothermal | 800 MW | ₱1,500 to ₱1,575 | 0.4 |
| Hydro | 2,500 MW | ₱300 | none |
| Oil-based | 1,500 MW | ₱14,000 to ₱14,700 | 3.0 |
| Solar | 2,500 MW | ₱40 | none |
| Wind | 400 MW | ₱60 | none |
 
**Constraints.**
 
| Constraint | Rows | Rule |
| --- | --- | --- |
| Demand | 6 | Total output in each block meets that block's demand |
| Capacity | 42 | Each group stays under its limit in each block. Solar is capped at 0, 250, 1,500, 1,375, 200 and 0 MW in B1 to B6. Wind is capped at 140, 120, 100, 100, 120 and 140 MW |
| Coal minimum | 6 | Coal runs at no less than 2,800 MW in every block |
| Coal ramping | 10 | Coal output changes by at most 1,500 MW between consecutive blocks, up and down |
| Hydro water | 1 | 4 × the sum of hydro output over all blocks is at most 8,000, a daily water budget |
| CO₂ cap | 1 | Total emissions stay under the limit you enter. Switched off unless you tick the box |
 
Unselected plant groups have their columns removed from the tableau. Any constraint that ends up with all-zero coefficients after that is dropped too, which is why Test Case #1 reports 36 unknowns and 59 constraints instead of 42 and 66.
 
## How it solves
 
1. `read_primal_csv()` loads the CSV and checks that it is numeric, that the last row is `Objective` and that the last column is `RHS`.
2. `build_initial_tableau()` keeps the selected groups, applies or drops the CO₂ row, then builds the dual maximization tableau `[ Aᵀ | I | 0 | c ]` over a bottom row of `[ -b | 0 | 1 | 0 ]`.
3. `Simplex()` pivots on the most negative entry of the bottom row and picks the pivot row by the smallest positive test ratio, using Gauss-Jordan elimination. It stops when the bottom row has no negative entries. If a pivot column has no positive entry, the problem is reported as infeasible.
4. `solve_func()` reads the dispatch x<sub>gt</sub> from the slack columns of the bottom row and the shadow prices S<sub>t</sub> from the Solution column. The price of electricity in block t is S<sub>t</sub> divided by the 4-hour block length.
The Simplex function can also return the tableau of every iteration, which is what feeds the Simplex Tableau tab.
 
## Test cases
 
| Case | Groups | CO₂ cap | Result |
| --- | --- | --- | --- |
| Test #1 | Coal, Natural Gas, Hydro, Oil-based, Solar, Wind | off | Optimal in 34 iterations, Z = ₱1,149,843,800.00. Prices in ₱/MWh for B1 to B6: 3,800.00, 6,565.00, 6,695.00, 14,700.00, 14,700.00, 14,280.00 |
| Test #2 | Coal, Natural Gas, Geothermal, Hydro, Oil-based, Wind | 180,000 | Optimal in 39 iterations, Z = ₱1,181,240,745.45. Prices in ₱/MWh for B1 to B6: 8,473.45, 8,538.45, 18,120.23, 18,400.23, 18,400.23, 8,603.45 |
| Test #3 | Coal, Natural Gas, Geothermal, Hydro, Oil-based, Wind | 165,000 | Infeasible after 37 iterations |
 
Test #3 uses the same plant groups as Test #2 with the cap lowered from 180,000 to 165,000 t, and that change alone turns the problem from solvable to infeasible.
 
## Input handling
 
- Pressing SOLVE with no plant group selected shows an "Invalid input" message.
- A CO₂ cap that is empty, not a number, or 0 or below shows an "Invalid input" message.
- If `data/Primal.csv` is missing or malformed, the app shows an error card instead of crashing.
## Files
 
| File | Purpose |
| --- | --- |
| `CerdeñaProject.R` | Simplex solver, tableau builder and the Shiny UI and server |
| `data/Primal.csv` | Primal tableau for all seven groups, 66 constraints and the objective row |
| `ui/styles.css` | Stylesheet for the whole interface |
| `ui/header.html` | Page header |
| `ui/custom_table.html` | Template shared by every table in the app |
| `docs/150-UserManual-Cerdeña.pdf` | User manual, English and Filipino |
 
