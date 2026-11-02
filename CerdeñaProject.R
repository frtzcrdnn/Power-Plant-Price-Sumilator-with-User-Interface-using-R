# Fritz Benedict C. Cerdeña | 2025 - 06390
# CMSC 150 Project - Luzon Grid Dispatch
#
# Pipeline:  primal tableau (CSV)  ->  pick groups / CO2 toggle  ->  transpose + slack + Z
#            ->  Simplex()  ->  Z, iterations, dispatch x_gt, prices S_t / h
#
# Run the test cases:   Rscript CerdeñaProject.R      (CSV must be in the working directory)


# CONSTANTS
CSV_PATH    <- "CMSC_150_-_Milestone_1_-_Primal_Tableau.csv" # primal tableau: ALL groups, CO2 cap ON (165,000)
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

# basic solution of a minimization tableau (bottom row), used for every iteration of the history
basic_solution_min <- function(tableau){
    numRows <- nrow(tableau); numCols <- ncol(tableau)
    bs <- tableau[numRows, 1:(numCols - 1)]
    bs[length(bs)] <- tableau[numRows, numCols] # the Z column shows Z
    return(bs)
}

# ------------------------------------------------------------------------------
# 5. ONE FUNCTION THAT DOES EVERYTHING (this is what the UI will call)
# ------------------------------------------------------------------------------
validate_input <- function(selected_groups, co2_on, co2_cap){
    if (length(selected_groups) == 0) return("No plant group selected.")
    if (isTRUE(co2_on) && (is.null(co2_cap) || length(co2_cap) != 1 || is.na(suppressWarnings(as.numeric(co2_cap))) || as.numeric(co2_cap) <= 0))
        return("The CO2 cap must be a number greater than 0.")
    return(NULL)
}

solve_dispatch <- function(selected_groups, co2_on = TRUE, co2_cap = 165000, primal = NULL, keepHistory = FALSE){
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
    out
}

# ------------------------------------------------------------------------------
# 6. HAND-OFF TO NEXT.JS: turn the result into JSON  (needs the jsonlite package:
#    install.packages("jsonlite"); it is only used to write JSON, not to optimize anything)
#    digits = NA keeps full precision (the default of 4 decimals would damage Z and the prices)
# ------------------------------------------------------------------------------
result_to_json <- function(res, include_history = FALSE){
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite: install.packages(\"jsonlite\")")
    keep <- c("status", "message", "groups", "co2_on", "co2_cap", "num_unknowns", "num_constraints", "iterations", "Z", "check_sum_bS")
    payload <- res[intersect(keep, names(res))]
    if (!is.null(res$dispatch)){
        payload$blocks   <- colnames(res$dispatch)
        payload$dispatch <- lapply(rownames(res$dispatch), function(g) as.list(res$dispatch[g, ]))
        names(payload$dispatch) <- rownames(res$dispatch)         # {"Coal": {"B1": 6160, ...}, ...}
        payload$prices   <- as.list(res$prices)                   # {"B1": 3800, ...}
    }
    if (include_history && !is.null(res$history)){
        payload$tableau_row_names <- rownames(res$initial_tableau)
        payload$tableau_col_names <- colnames(res$initial_tableau)
        payload$history <- lapply(res$history, function(h) list(
            iteration = h$iteration,
            pivotRow  = if (is.na(h$pivotRow)) NA_integer_ else as.integer(h$pivotRow),   # 1-based; null when there is no pivot
            pivotCol  = if (is.na(h$pivotCol)) NA_integer_ else as.integer(h$pivotCol),
            tableau   = unname(h$tableau),                              # array of rows
            basic_solution = as.list(basic_solution_min(h$tableau))))
    }
    jsonlite::toJSON(payload, auto_unbox = TRUE, digits = NA, na = "null", pretty = FALSE)
}

# ------------------------------------------------------------------------------
# 7. TEST CASES from the Project Specs
# ------------------------------------------------------------------------------
print_result <- function(name, res){
    cat("=====", name, "=====\n")
    cat("Plant groups :", paste(res$groups, collapse = ", "), "\n")
    cat("CO2 cap      :", if (isTRUE(res$co2_on)) paste0("ON (", format(res$co2_cap, big.mark = ","), " t/day)") else "OFF", "\n")
    if (res$status == "invalid"){ cat("Invalid input:", res$message, "\n\n"); return(invisible()) }
    cat("LP setup     :", res$num_unknowns, "unknowns,", res$num_constraints, "constraints\n")
    cat("Iterations   :", res$iterations, "\n")
    if (res$status == "infeasible"){ cat("Result       : The problem is infeasible.\n\n"); return(invisible()) }
    cat("Total minimum cost, Z :", format(res$Z, big.mark = ",", nsmall = 2), "\n")
    cat("Unknowns x_gt (MW):\n");              print(round(res$dispatch, 2))
    cat("Price = S_t / 4 (PHP/MWh):\n");       print(round(res$prices, 2))
    cat("Check sum(b * S) - Z =", format(res$check_sum_bS - res$Z, digits = 3), "\n\n")
}

check_case <- function(case, res){
    ok <- res$status == case$status && res$iterations == case$iterations
    if (case$status == "optimal"){
        ok <- ok && abs(res$Z - case$Z) / case$Z < 1e-4 &&                      # within 0.01 %
                    all(abs(res$prices - case$prices) <= 0.01) &&               # within 0.01 PHP/MWh
                    abs(res$check_sum_bS - res$Z) / res$Z < 1e-9
    }
    cat(if (ok) "PASS" else "FAIL", "-", case$name, "\n\n")
    ok
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

if (sys.nframe() == 0){ # only runs with `Rscript CerdeñaProject_CSV.R`, not when this file is source()d
    primal <- read_primal_csv(CSV_PATH)
    cat("CSV read:", nrow(primal$A), "constraints x", ncol(primal$A), "unknowns\n\n")
    passed <- 0
    for (case in TEST_CASES){
        res <- solve_dispatch(case$groups, case$co2_on, case$co2_cap, primal)
        print_result(case$name, res)
        if (check_case(case, res)) passed <- passed + 1
    }
    cat(passed, "of", length(TEST_CASES), "test cases passed\n")

    cat("\n--- Invalid input examples ---\n")
    print_result("No group selected", c(solve_dispatch(character(0), TRUE, 165000, primal), list(groups = "none")))
    print_result("CO2 cap of 0",      c(solve_dispatch("Coal", TRUE, 0, primal), list(groups = "Coal", co2_on = TRUE, co2_cap = 0)))
}
