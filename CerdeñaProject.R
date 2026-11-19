# Fritz Benedict C. Cerdeña | 2025 - 06390

# Pipeline:  primal tableau (in CSV)  ->  pick groups / CO2 toggle  ->  transpose + slack + Z
# ->  Simplex()  ->  Z, iterations, dispatch x_gt, prices S_t / h   ->  R Shiny UI

# NEEDS: 
# 1) install.packages("shiny")
# 2) put "data/Primal.csv" ->, "ui/styles.css", "ui/header.html" and "ui/custom_table.html" in their respective folder

library(shiny)

# CONSTANTS
CSV_PATH    <- "data/Primal.csv" # primal tableau: ALL groups, CO2 cap ON (165,000)
BLOCK_HOURS <- 4                                              # h: hours per block (prices = S_t / h)
N_BLOCKS    <- 6

# UI / display name -> prefix used in the CSV column names (Coal_B1, NatGas_B1, Geo_B1, ...)
GROUP_KEYS <- c("Coal" = "Coal", "Natural Gas" = "NatGas", "Geothermal" = "Geo", "Hydro" = "Hydro",
                "Oil-based" = "Oil", "Solar" = "Solar", "Wind" = "Wind")

# keepHistory = TRUE also returns $history: the tableau of every iteration, defaults to FALSE
Simplex <- function(tableau, isMax, keepHistory = FALSE){
  numRows <- nrow(tableau)
  numCols <- ncol(tableau)
  tol <- 1e-8 # for handling near 0 values
  
  iterations <- 0 # count pivots
  # ADDED: For project specifications, infeasible is needed
  infeasible <- FALSE 
  history <- list() # NEW: tableau of every iteration (only filled when keepHistory = TRUE)
  while(TRUE){
    bottomRow <- tableau[numRows, 1:(numCols - 1)] # gets the bottom row, from first column to the 2nd last column
    
    # checks if there are no more negative numbers left; end the loop
    if (!any(bottomRow < -tol)){
      if (keepHistory) history[[length(history) + 1]] <- list(iteration = iterations, pivotRow = NA, pivotCol = NA, tableau = tableau) # NEW: final tableau
      break
    }
    
    minVal <- min(bottomRow)
    pivotCol <- which(abs(bottomRow - minVal) < tol)[1] # returns which column is equal to the minimum value
    
    validRows <- c()    
    testRatios <- c()
    
    # get all the test ratios
    for (i in 1:(numRows - 1)){ # loops through all rows except for bottom row
      b <- tableau[i, pivotCol]
      if (b > tol){ # checks if positive, and not dividing by 0
        a <- tableau[i, numCols]
        tr <- a/b
        validRows <- c(validRows, i) # stores the row index if the row is valid as pivot row
        testRatios <- c(testRatios, tr) # stores the test ratios for later
      }
    }
    
    # if no valid pivot rows
    if (length(validRows) == 0){
      infeasible <- TRUE
      if (keepHistory) history[[length(history) + 1]] <- list(iteration = iterations, pivotRow = NA, pivotCol = pivotCol, tableau = tableau) # NEW
      break
    }
    
    # get the smallest test ratio
    minTr <- min(testRatios)
    pivotRow <- validRows[which(abs(testRatios - minTr) < tol)][1] # checks which of the testratios would result in a difference of 0
    PivotElement <- tableau[pivotRow, pivotCol]
    
    if (keepHistory) history[[length(history) + 1]] <- list(iteration = iterations, pivotRow = pivotRow, pivotCol = pivotCol, tableau = tableau) # NEW: state BEFORE this pivot
    iterations <- iterations + 1
    # (gauss-jordan / row normalization) and then elimination/update
    tableau[pivotRow, ] <- tableau[pivotRow, ] / PivotElement # nPR = PR/PE
    for (i in 1:numRows){
      if (i != pivotRow){
        C <- tableau[i, pivotCol] # value to eliminate
        tableau[i, ] <- tableau[i, ] - (tableau[pivotRow, ] * C) # Rn = Rn - (nPR * C)
      }
    }
  }
  
  finalTableau <- tableau
  Z <- unname(tableau[numRows, numCols]) # unname() removes the col/row names
  colNames <- colnames(tableau) # gets the column names from tableau
  variables <- colNames[1:(numCols - 1)] # -1 since sobra ng isa
  
  # if we are finding the maximum
  if(isMax){
    basicSolution <- numeric(numCols - 1)
    names(basicSolution) <- variables
    for (i in 1:(numCols - 1)){ # loops through all columns except the solution column
      colVals <- tableau[, i]
      nonZeroRows <- which(abs(colVals) > tol) # saves rows indexes in which the value is not 0
      if (length(nonZeroRows) == 1 && abs(colVals[nonZeroRows] - 1) < tol){ # verifies if there is only 1 nonzero per column, and ensures that it is 1
        basicSolution[i] <- tableau[nonZeroRows, numCols] / tableau[nonZeroRows, i] # assigns it to the basic solutions
      }else{
        basicSolution[i] <- 0 # if there are more than 1 nonzero in the column, assign 0 to the basic solution
      }
    }
  }else{ # if finding minimum
    # the basic solution is already right at the bottom row
    basicSolution <- tableau[numRows, 1:(numCols - 1)]
    names(basicSolution) <- variables
    basicSolution[length(basicSolution)] <- Z
  }
  
  return(list(finalTableau = finalTableau, basicSolution = basicSolution, Z = Z, iterations = iterations, infeasible = infeasible, history = history))
}


# Reads the Primal Tableau from CSV File to be used in building the initial tableau
read_primal_csv <- function(path = CSV_PATH){
  if (!file.exists(path)) stop("Cannot find the primal tableau CSV: ", path)
  
  # reads the .csv file and converts to M matrix in R
  # https://stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html
  df <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  M  <- as.matrix(df)
  
  # input validations
  if (!is.numeric(M) || anyNA(M)) stop("The CSV must contain only numbers (besides row/column names).")
  if (tolower(rownames(M)[nrow(M)]) != "objective") stop("The last row of the CSV must be the Objective row.")
  if (tolower(colnames(M)[ncol(M)]) != "rhs") stop("The last column of the CSV must be RHS.")
  
  var_names <- colnames(M)[-ncol(M)]
  # https://shiny.abdn.ac.uk/Stats/tutorials/R_Intermediate/R_Intermediate_P1/matrix-subsetting.html
  list(A = M[-nrow(M), -ncol(M), drop = FALSE],   # constraint coefficients
       b = M[-nrow(M), ncol(M)],                  # RHS of the constraints
       c = M[nrow(M), -ncol(M)],                  # objective row (h * c_gt)
       var_names = var_names,
       row_names = rownames(M)[-nrow(M)],
       groups    = unique(sub("_B[0-9]+$", "", var_names))) # extract sub-groups (unique), (Coal_B1 etc. --> Coal)
}

# CSV prefix -> display name (what the UI and the tables show)
display_name <- function(keys) unname(ifelse(keys %in% GROUP_KEYS, names(GROUP_KEYS)[match(keys, GROUP_KEYS)], keys))

# accepts display names ("Natural Gas") or CSV prefixes ("NatGas")
resolve_groups <- function(selected_groups, available){
  keys <- ifelse(selected_groups %in% names(GROUP_KEYS), GROUP_KEYS[selected_groups], selected_groups) # specifically used ifelse for vectors
  unknown <- setdiff(keys, available) # var_names present in keys but not in available are assigned to unknown var_names
  if (length(unknown) > 0) stop("Unknown plant group: ", paste(unknown, collapse = ", "))
  return(unname(keys))
}


# Builds the Initial Tableau
# - remove unselected columns
# - toggle setting for CO2_Cap ON/OFF
# - remove rows that are now zeroes
build_initial_tableau <- function(primal, selected_groups, co2_on = TRUE, co2_cap = 165000){
  keys <- resolve_groups(selected_groups, primal$groups)
  
  # columns of the selected groups only (keeps the CSV order, not the order the user clicked)
  # RegEx used: https://search.r-project.org/R/refmans/base/html/regex.html
  keep_cols <- sub("_B[0-9]+$", "", primal$var_names) %in% keys # sub replaces all matches 
  A         <- primal$A[, keep_cols, drop = FALSE]
  C_primal  <- primal$c[keep_cols]
  var_names <- primal$var_names[keep_cols]
  b_vals    <- primal$b
  row_names <- primal$row_names
  
  # if CO2 cap OFF -> drop the row, else ON -> use the limit typed by the user
  is_co2 <- row_names == "CO2_Cap"
  keep_rows <- rep(TRUE, length(row_names))
  if (co2_on) b_vals[is_co2] <- -co2_cap else keep_rows[is_co2] <- FALSE
  
  # drop constraints that belong only to removed groups that are all zeros
  keep_rows <- keep_rows & (rowSums(abs(A)) > 0)
  A         <- A[keep_rows, , drop = FALSE]
  b_vals    <- unname(b_vals[keep_rows])
  row_names <- row_names[keep_rows]
  
  num_vars        <- ncol(A)  # n
  num_constraints <- nrow(A)  # m
  
  # [ A^T | I_n | 0_z | C_primal ]
  top_block  <- cbind(t(A), diag(num_vars), matrix(0, num_vars, 1), matrix(C_primal, num_vars, 1))
  # [ -b_primal | 0_slack | 1_z | 0_sol ]
  bottom_row <- c(-b_vals, numeric(num_vars), 1, 0)
  
  initial_tableau <- rbind(top_block, bottom_row)
  
  # descriptive names: S1_Demand_B1 ... | Coal_B1 ... (slack = the original unknown) | Z | Solution
  colnames(initial_tableau) <- c(paste0("S", 1:num_constraints, "_", row_names), var_names, "Z", "Solution")
  rownames(initial_tableau) <- c(var_names, "Z")
  
  return(list(tableau = initial_tableau, m = num_constraints, n = num_vars, b = b_vals, constraint_names = row_names, var_names = var_names, groups = keys))
}

# shadow prices from final tableau
shadow_prices <- function(final_tableau, m, tol = 1e-8){
  numRows <- nrow(final_tableau)
  numCols <- ncol(final_tableau)
  S <- numeric(m)
  for (i in 1:m){
    col <- final_tableau[, i]
    nz  <- which(abs(col) > tol) # not zero
    if (length(nz) == 1 && nz < numRows && abs(col[nz] - 1) < tol) S[i] <- final_tableau[nz, numCols]
  }
  names(S) <- colnames(final_tableau)[1:m]
  return(S)
}

# basic solution of a minimization tableau (bottom row), used for every iteration
basic_solution_min <- function(tableau){
  numRows <- nrow(tableau); numCols <- ncol(tableau)
  bs <- tableau[numRows, 1:(numCols - 1)]
  bs[length(bs)] <- tableau[numRows, numCols] # the Z column shows Z
  return(bs)
}

# func for input validations
validate_input <- function(selected_groups, co2_on, co2_cap){
  if (length(selected_groups) == 0) return("No plant group selected.")
  if (isTRUE(co2_on) && (is.null(co2_cap) || length(co2_cap) != 1 || is.na(suppressWarnings(as.numeric(co2_cap))) || as.numeric(co2_cap) <= 0))
    return("The CO2 cap must be a number greater than 0.")
  return(NULL)
}

solve_func <- function(selected_groups, co2_on = TRUE, co2_cap = 165000, primal = NULL, keepHistory = FALSE){
  msg <- validate_input(selected_groups, co2_on, co2_cap)
  if (!is.null(msg)) return(list(status = "invalid", message = msg))
  
  if (is.null(primal)) primal <- read_primal_csv(CSV_PATH)
  if (co2_on) co2_cap <- as.numeric(co2_cap)
  built <- build_initial_tableau(primal, selected_groups, co2_on, co2_cap)
  res   <- Simplex(built$tableau, FALSE, keepHistory)
  
  out <- list(status = NA, message = NA, groups = display_name(built$groups), co2_on = co2_on, co2_cap = if (co2_on) co2_cap else NA,
              num_unknowns = built$n, num_constraints = built$m, iterations = res$iterations,
              initial_tableau = built$tableau)
  
  if (res$infeasible){
    out$status  <- "infeasible"
    out$message <- "The problem is infeasible."
  } else {
    # dispatch x_gt = slack columns of the BOTTOM ROW; prices = S_t / h from the Solution column
    x <- res$basicSolution[built$m + 1:built$n]
    x[abs(x) < 1e-9] <- 0
    dispatch <- matrix(x, nrow = length(built$groups), ncol = N_BLOCKS, byrow = TRUE,
                       dimnames = list(display_name(built$groups), paste0("B", 1:N_BLOCKS)))
    S <- shadow_prices(res$finalTableau, built$m)
    
    out$status   <- "optimal"
    out$message  <- "Optimal solution found"
    out$Z        <- res$Z
    out$dispatch <- dispatch                                   # rows = groups, cols = B1..B6 (MW)
    out$prices   <- setNames(S[1:N_BLOCKS] / BLOCK_HOURS, paste0("B", 1:N_BLOCKS)) # PHP/MWh
    out$check_sum_bS <- sum(built$b * S)                       # must equal Z
    out$final_tableau <- res$finalTableau
  }

  if (keepHistory) out$history <- res$history
  return(out)
}

TEST_CASES <- list(
  list(name = "Test Case #1", groups = c("Coal", "Natural Gas", "Hydro", "Oil-based", "Solar", "Wind"),
       co2_on = FALSE, co2_cap = 165000, status = "optimal", iterations = 34, Z = 1149843800,
       prices = c(3800, 6565, 6695, 14700, 14700, 14280)),
  list(name = "Test Case #2", groups = c("Coal", "Natural Gas", "Geothermal", "Hydro", "Oil-based", "Wind"),
       co2_on = TRUE, co2_cap = 180000, status = "optimal", iterations = 39, Z = 1181240745.45,
       prices = c(8473.45, 8538.45, 18120.23, 18400.23, 18400.23, 8603.45)),
  list(name = "Test Case #3", groups = c("Coal", "Natural Gas", "Geothermal", "Hydro", "Oil-based", "Wind"),
       co2_on = TRUE, co2_cap = 165000, status = "infeasible", iterations = 37)
)

# =====================================================================================
#                                    R SHINY USER INTERFACE
# =====================================================================================

PESO <- "\u20b1"  # peso sign 
# display-only info for the UI (the numbers used by the solver all come from Primal.csv)
GROUP_COLORS <- c("Coal" = "#4b5563", "Natural Gas" = "#f08a24", "Geothermal" = "#a66bbe", "Hydro" = "#2f6fb3",
                  "Oil-based" = "#8b5a3c", "Solar" = "#f2c200", "Wind" = "#3fa66a")
GROUP_CAP    <- c("Coal" = 7000, "Natural Gas" = 3300, "Geothermal" = 800, "Hydro" = 2500,
                  "Oil-based" = 1500, "Solar" = 2500, "Wind" = 400)
BLOCK_LABELS <- c("12\u20134 AM", "4\u20138 AM", "8 AM\u201312 NN", "12\u20134 PM", "4\u20138 PM", "8 PM\u201312 MN")

# NUMBER FORMATTING
# trims trailing zeros: 7000.00 -> "7,000" ; 4540.909 -> "4,540.91"
fmt_num <- function(x, digits = 2){
  x   <- ifelse(abs(x) < 1e-9, 0, x)
  out <- formatC(x, format = "f", digits = digits, big.mark = ",") # big.mark makes 123456.789 -> 123,456.79
  out[out %in% c("-0", "")] <- "0" # rounding sometimes cause -0, turns to 0
  return(out)
}
# always 2 decimals: 3800 -> "3,800.00"
fmt_fixed <- function(x) formatC(ifelse(abs(x) < 1e-9, 0, x), format = "f", digits = 2, big.mark = ",")

esc <- function(x) htmltools::htmlEscape(x)

# ---- small HTML builders (plain tables so we control scrolling, sticky headers and highlights) ----
alert_box <- function(type, title, body = NULL){
  div(class = paste("alert-box", type), div(class = "alert-title", title), if (!is.null(body)) div(class = "alert-body", body))
}

metric_card <- function(label, value, sub = NULL, cls = ""){
  div(class = paste("metric", cls), div(class = "metric-label", label), div(class = "metric-value", value),
      if (!is.null(sub)) div(class = "metric-sub", sub))
}

# shared table wrapper: the markup lives in custom_table.html, R only fills the placeholders
# style = inline CSS for the scroll box, thead / tbody = the <tr> rows (HTML strings)
table_component <- function(thead, tbody, maxh = NULL){
  htmlTemplate("ui/custom_table.html",
               style = HTML(if (is.null(maxh)) "" else paste0("max-height:", maxh, "px")),
               thead = HTML(thead), tbody = HTML(tbody),
               document_ = FALSE)
}

# generic scrollable matrix -> HTML. cls = matrix of css classes for each cell, hcls = css class of each column header
grid_html <- function(M, rn, cn, cls = NULL, hcls = NULL, corner = "", maxh = 520){
  nr <- nrow(M); nc <- ncol(M)
  if (is.null(cls))  cls  <- matrix("", nr, nc)
  if (is.null(hcls)) hcls <- rep("", nc)
  head <- paste0("<tr><th class='corner'>", esc(corner), "</th>",
                 paste0("<th class='", hcls, "'>", esc(cn), "</th>", collapse = ""), "</tr>")
  body <- vapply(seq_len(nr), function(i){
    paste0("<tr><th class='rh'>", esc(rn[i]), "</th>",
           paste0("<td class='", cls[i, ], "'>", M[i, ], "</td>", collapse = ""), "</tr>")
  }, character(1))
  table_component(head, paste(body, collapse = ""), maxh)
}

# Dispatch table x_gt (MW) with a colour-intensity heat shading per row
dispatch_html <- function(d){
  nr <- nrow(d)
  cells <- vapply(seq_len(nr), function(i){
    g    <- rownames(d)[i]
    rgbv <- col2rgb(GROUP_COLORS[g])
    mx   <- max(d[i, ])
    a    <- if (mx > 0) 0.06 + 0.44 * d[i, ] / mx else rep(0, ncol(d))
    paste0("<tr><th class='rh'><span class='dot' style='background:", GROUP_COLORS[g], "'></span>", esc(g), "</th>",
           paste0("<td style='background:rgba(", rgbv[1], ",", rgbv[2], ",", rgbv[3], ",", sprintf("%.2f", a), ")'>",
                  fmt_num(d[i, ]), "</td>", collapse = ""), "</tr>")
  }, character(1))
  total <- paste0("<tr class='total'><th class='rh'>Total</th>", paste0("<td>", fmt_num(colSums(d)), "</td>", collapse = ""), "</tr>")
  head  <- paste0("<tr><th class='corner'>Plant group</th>",
                  paste0("<th>", colnames(d), "<small>", BLOCK_LABELS, "</small></th>", collapse = ""), "</tr>")
  table_component(head, paste0(paste(cells, collapse = ""), total))
}

prices_html <- function(prices, h = BLOCK_HOURS){
  head <- paste0("<tr><th class='corner'></th>", paste0("<th>", names(prices), "<small>", BLOCK_LABELS, "</small></th>", collapse = ""), "</tr>")
  r1   <- paste0("<tr><th class='rh'>Shadow price S<sub>t</sub> (", PESO, ")</th>", paste0("<td>", fmt_fixed(prices * h), "</td>", collapse = ""), "</tr>")
  r2   <- paste0("<tr class='price-row'><th class='rh'>Price = S<sub>t</sub> / ", h, " (", PESO, "/MWh)</th>",
                 paste0("<td>", fmt_fixed(prices), "</td>", collapse = ""), "</tr>")
  table_component(head, paste0(r1, r2))
}

# non-zero shadow prices of ALL constraints
shadow_html <- function(res){
  m   <- res$num_constraints
  S   <- shadow_prices(res$final_tableau, m)
  b   <- -unname(res$initial_tableau[nrow(res$initial_tableau), 1:m])   # bottom row of the initial tableau holds -b
  nm  <- sub("^S[0-9]+_", "", names(S))
  idx <- which(abs(S) > 1e-9)
  if (length(idx) == 0) return(HTML("<p class='muted'>All shadow prices are 0.</p>"))
  rows <- vapply(idx, function(i) paste0("<tr><th class='rh'>S<sub>", i, "</sub> &nbsp;", esc(nm[i]), "</th><td>", fmt_fixed(S[i]),
                                         "</td><td>", fmt_num(b[i]), "</td><td>", fmt_fixed(b[i] * S[i]), "</td></tr>"), character(1))
  head <- paste0("<tr><th class='corner'>Constraint (non-zero S<sub>i</sub> only)</th>",
                 "<th>S<sub>i</sub></th><th>RHS b<sub>i</sub></th><th>b<sub>i</sub> &middot; S<sub>i</sub></th></tr>")
  table_component(head, paste0(paste(rows, collapse = "")), maxh = 320)
}

# ---- UI ----------------------------------------------------------------------------

ui <- fluidPage(
  title = "Luzon Grid Dispatch Simulator",
  tags$head(includeCSS("ui/styles.css")),   # stylesheet: styles.css
  includeHTML("ui/header.html"),            # header component: header.html
  div(class = "app-body",
      sidebarLayout(
        sidebarPanel(width = 4,
                     div(class = "card-lite",
                         h4(span(class = "step-no", "1"), "Plant groups"),
                         checkboxGroupInput("groups", NULL,
                                            choiceNames  = lapply(names(GROUP_KEYS), function(g)
                                              HTML(paste0("<span class='dot' style='background:", GROUP_COLORS[g], "'></span>", g,
                                                          "<span class='cap'>", format(GROUP_CAP[g], big.mark = ","), " MW</span>"))),
                                            choiceValues = names(GROUP_KEYS)),
                         div(class = "btn-row",
                             actionButton("check_all", "Check All", class = "btn-soft"),
                             actionButton("reset", "Reset", class = "btn-soft"))),
                     div(class = "card-lite",
                         h4(span(class = "step-no", "2"), "CO\u2082 cap"),
                         checkboxInput("co2_on", "Limit daily CO\u2082 emissions", value = FALSE),
                         conditionalPanel("input.co2_on",
                                          numericInput("co2_cap", "Limit (t CO\u2082 / day)", value = 165000, min = 0, step = 1000))),
                     div(class = "card-lite",
                         h4(span(class = "step-no", "3"), "Solve"),
                         actionButton("solve", HTML("SOLVE")),
                         hr(style = "margin:14px 0 8px"),
                         p(class = "muted", style = "font-size:12.5px; margin:0 0 6px", "Pre-loaded Test Cases:"),
                         div(class = "btn-row test-row",
                             actionButton("tc1", "Test #1", class = "btn-soft"),
                             actionButton("tc2", "Test #2", class = "btn-soft"),
                             actionButton("tc3", "Test #3", class = "btn-soft")))
        ),
        mainPanel(width = 8,
                  tabsetPanel(id = "tabs",
                              tabPanel("Results", uiOutput("results_ui")),
                              tabPanel("Simplex Tableau", uiOutput("tableau_tab"))
                  )
        )
      )
  )
)

# ---- Server ------------------------------------------------------------------------
server <- function(input, output, session){
  # the primal tableau is read ONCE; any problem is turned into a friendly message (never a raw error)
  primal <- tryCatch(read_primal_csv(CSV_PATH), error = function(e) NULL)
  rv <- reactiveValues(res = NULL)
  
  # ---- inputs: Check All / Reset / sample test cases ----
  observeEvent(input$check_all, updateCheckboxGroupInput(session, "groups", selected = names(GROUP_KEYS)))
  observeEvent(input$reset, {
    updateCheckboxGroupInput(session, "groups", selected = character(0))
    updateCheckboxInput(session, "co2_on", value = FALSE)
    updateNumericInput(session, "co2_cap", value = 165000)
    rv$res <- NULL
  })
  load_case <- function(idx){
    tc <- TEST_CASES[[idx]]
    updateCheckboxGroupInput(session, "groups", selected = tc$groups)
    updateCheckboxInput(session, "co2_on", value = tc$co2_on)
    updateNumericInput(session, "co2_cap", value = tc$co2_cap)
  }
  observeEvent(input$tc1, load_case(1))
  observeEvent(input$tc2, load_case(2))
  observeEvent(input$tc3, load_case(3))
  
  # ---- solve ----
  observeEvent(input$solve, {
    if (is.null(primal)){
      rv$res <- list(status = "error", message = paste0("The primal tableau file \"", CSV_PATH, "\" could not be read. Place it in the same folder as the app and restart."))
    } else {
      rv$res <- tryCatch(
        solve_func(input$groups %||% character(0), isTRUE(input$co2_on), input$co2_cap, primal, keepHistory = TRUE),
        error = function(e) list(status = "error", message = "Something went wrong while solving. Please check your selection and try again."))
    }
    updateTabsetPanel(session, "tabs", selected = "Results")
  })
  
  # ---- Results tab ----
  output$results_ui <- renderUI({
    res <- rv$res
    if (is.null(res)){
      return(div(class = "card-lite welcome",
                 h4("Welcome!"),
                 tags$ol(
                   tags$li("Tick the plant groups that are available (or press Check All)."),
                   tags$li("Switch on the CO\u2082 cap and edit its limit. (Optional)"),
                   tags$li("Press ", tags$b("SOLVE"), ".")),
                 p(class = "muted", "You will get the minimum daily cost, the output of every plant in every 4-hour block, and the price of electricity in each block. ",
                   "The ", tags$b("Simplex Tableau"), " tab lets you see every iteration.")))
    }
    if (res$status == "invalid") return(alert_box("warn", "Invalid input", res$message))
    if (res$status == "error")   return(alert_box("bad", "Unable to solve", res$message))
    if (res$status == "infeasible"){
      return(tagList(
        alert_box("bad", "The problem is infeasible.",
                  paste0("The selected plants cannot meet demand in every block under these limits (no test ratio could be chosen). Number of iterations: ", res$iterations, ".")),
        div(class = "metrics", metric_card("Iterations", res$iterations), metric_card("Unknowns", res$num_unknowns), metric_card("Constraints", res$num_constraints)),
        p(class = "muted", "Tip: select more plant groups, or raise / switch off the CO\u2082 cap. The Simplex Tableau tab shows where the algorithm got stuck.")))
    }
    
    tagList(
      alert_box("ok", "Optimal solution found"),
      div(class = "metrics",
          metric_card("Total minimum cost, Z", paste0(PESO, fmt_fixed(res$Z)), "for one day", "main"),
          metric_card("Iterations", res$iterations),
          metric_card("Unknowns", res$num_unknowns), # Xnt variables
          metric_card("Constraints", res$num_constraints)),
      div(class = "card-lite", h4("Output of each plant group"), dispatch_html(res$dispatch)),
      div(class = "card-lite", h4("Price of electricity"), prices_html(res$prices)),
      div(class = "card-lite",
          h4("Shadow prices & check"),
          shadow_html(res))
            )
  })
  
  demand <- reactive(if (is.null(primal)) NULL else unname(primal$b[paste0("Demand_B", 1:N_BLOCKS)]))
  
  # ---- Simplex Tableau tab ----
  output$tableau_tab <- renderUI({
    res <- rv$res
    if (is.null(res) || is.null(res$history)){
      return(alert_box("bad", "No tableau yet", "Press SOLVE to get the iteratons of the Simplex Tableau."))
    }
    N <- length(res$history) - 1
    tagList(
      div(class = "card-lite",
          div(class = "iter-ctrl",
              actionButton("iter_first", HTML("&#9198; First"), class = "btn-soft"),
              actionButton("iter_prev",  HTML("&#9664; Prev"),  class = "btn-soft"),
              actionButton("iter_next",  HTML("Next &#9654;"),  class = "btn-soft"),
              actionButton("iter_last",  HTML("Last &#9197;"),  class = "btn-soft")),
          sliderInput("iter", "Iteration", min = 0, max = N, value = 0, step = 1, width = "100%",
                      animate = animationOptions(interval = 1000, loop = FALSE))),
      uiOutput("tableau_view"))
  })
  
  # iteration slider
  set_iter <- function(v) updateSliderInput(session, "iter", value = v)
  observeEvent(input$iter_first, set_iter(0))
  observeEvent(input$iter_prev,  set_iter(max(0, (input$iter %||% 0) - 1)))
  observeEvent(input$iter_next,  { req(rv$res$history); set_iter(min(length(rv$res$history) - 1, (input$iter %||% 0) + 1)) })
  observeEvent(input$iter_last,  { req(rv$res$history); set_iter(length(rv$res$history) - 1) })
  
  # initial matrices
  output$tableau_view <- renderUI({
    res <- rv$res; req(!is.null(res), !is.null(res$history))
    hist <- res$history
    k    <- min(max((input$iter %||% 0), 0), length(hist) - 1) + 1
    h    <- hist[[k]]; tab <- h$tableau
    nr <- nrow(tab); nc <- ncol(tab); m <- res$num_constraints; n <- res$num_unknowns
    pr <- h$pivotRow; pc <- h$pivotCol
    
    #   
    info <- if (!is.na(pr)){ 
      pe <- tab[pr, pc]
      paste0("Pivot column: <b>", esc(colnames(tab)[pc]), "</b> &nbsp;\u00b7&nbsp; ",
             "Pivot row: <b>", esc(rownames(tab)[pr]), "</b> &nbsp;\u00b7&nbsp; ",
             "Pivot element: <b>", fmt_num(pe, 4), "</b>")
    } else if (!is.na(pc)){
      paste0("<b>Infeasible.</b> Pivot column <b>", esc(colnames(tab)[pc]), "</b> has no positive entry, so there is no test ratio to choose from.")
    } else {
      "<b>Final tableau.</b> No negative entries are left in the bottom row, so the solution is optimal."
    }
    
    # cell classes: pivot row / column / cell, negative bottom-row entries, zeros
    cls <- matrix("", nr, nc)
    cls[abs(tab) < 1e-9] <- "zero"
    cls[nr, which(tab[nr, 1:(nc - 1)] < -1e-8)] <- "neg"
    if (!is.na(pr)) cls[pr, ] <- "prow"
    if (!is.na(pc)) cls[, pc] <- ifelse(cls[, pc] == "neg", "neg pcol", "pcol")
    if (!is.na(pr) && !is.na(pc)) cls[pr, pc] <- "pcell"
    hcls <- c(rep("hS", m), rep("hX", n), "hZ", "hZ")
    
    # basic solution of this iteration: the bottom row, with the Z column showing Z
    bs <- basic_solution_min(tab)
    bs_cls <- matrix(ifelse(abs(bs) > 1e-9, "nz", "zero"), nrow = 1)
    
    tagList(
      div(class = "pivot-info", HTML(paste0("<b>Iteration ", h$iteration, " of ", length(hist) - 1, "</b> &nbsp;|&nbsp; ", info))),
      div(class = "card-lite",
          h4("Simplex tableau"),
          grid_html(matrix(fmt_num(tab, 4), nr, nc), rownames(tab), colnames(tab), cls, hcls, corner = "Row \\ Column", maxh = 560),
          p(class = "legend-line",
            tags$span(class = "sw", style = "background:#7a3b0c"), "Shadow Prices",
            tags$span(class = "sw", style = "background:#0f2540"), "Slack/Unknowns",
            tags$span(class = "sw", style = "background:#fff3d6;border:1px solid #ddd"), "Pivot Row",
            tags$span(class = "sw", style = "background:#e1effd;border:1px solid #ddd"), "Pivot Column",
            tags$span(class = "sw", style = "background:#f5a623"), "Pivot Element")),
      div(class = "card-lite",
          h4("Basic solution"),
          grid_html(matrix(fmt_num(unname(bs), 4), 1, length(bs)), "Basic Solution", names(bs), bs_cls, hcls[-length(hcls)], corner = "", maxh = 140)))
  })
}

# =====================================================================================
#                                       LAUNCH PROGRAM
# =====================================================================================

# can launch directly using Ctrl + Shift + s
app <- shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))