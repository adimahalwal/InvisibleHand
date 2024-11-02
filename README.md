# InvisibleHand

Monte Carlo simulations 

learn more about it: https://youtu.be/r7cn3WS5x9c?si=ZuaV7qInI_g_yc_D

### Find out more about building applications with Shiny here:
###    https://shiny.posit.co/


### For building interactive web applications.
library(shiny)
### A collection of data manipulation and visualization tools.
library(tidyverse)
### A collection of data manipulation and visualization tools.
library(MASS)
### To deploy applications to RStudio's hosting platform.
library(rsconnect) 

## Define UI

**UI Definition: ui <- fluidPage(...)

**fluidPage() creates a responsive web page layout for the app.
Title: titlePanel("Monte Carlo Simulation with Different OLS Procedures")

**Displays the title at the top of the app.
Layout: sidebarLayout(...)

**Splits the UI into a sidebar for inputs and a main panel for outputs.
Sidebar Panel (sidebarPanel(...)):

**Simulation Settings (h4("Simulation Settings")): Headers and input elements for setting up simulations.

sliderInput("n", ...): Slider for selecting sample size, with a range from 10 to 500.
sliderInput("num_sim", ...): Slider for setting the number of simulations (100 to 5000).
sliderInput("sigma", ...): Slider to adjust error standard deviation (1 to 20).
Model Parameters (h4("Model Parameters")):

**selectInput("procedure", ...): Dropdown to select one of three OLS procedures.
Conditional Inputs (conditionalPanel(...)):
Only appears if a procedure is selected.
numericInput("beta_0", ...): Input for the intercept value.
numericInput("beta_1", ...) and numericInput("beta_2", ...): Inputs for slope values.
sliderInput("rho", ...): Slider for correlation, ranging from -1 to 1.
Main Panel (mainPanel(...)):

**Outputs:
verbatimTextOutput("results"): Text output for displaying simulation results.
plotOutput("density_plots"): Plot area to show density plots of results.


ui <- fluidPage(
  titlePanel("Monte Carlo Simulation with Different OLS Procedures"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Simulation Settings"),
      sliderInput("n", "Sample Size (n):", min = 10, max = 500, value = 50),
      sliderInput("num_sim", "Number of Simulations:", min = 100, max = 5000, value = 1000),
      sliderInput("sigma", "Error Standard Deviation (sigma):", min = 1, max = 20, value = 10),
      
      h4("Model Parameters"),
      selectInput("procedure", "Procedure:", choices = c("Procedure 1" = "procedure1", "Procedure 2" = "procedure2", "Procedure 3" = "procedure3")),
      
      conditionalPanel(
        condition = "input.procedure == 'procedure1' || input.procedure == 'procedure2'|| input.procedure == 'procedure3'",
        numericInput("beta_0", "True Intercept (beta_0):", 2),
        numericInput("beta_1", "True Slope (beta_1):", 1),
        numericInput("beta_2", "True Slope (beta_2):", 1),
        sliderInput("rho", "Correlation (rho):", min = -1, max = 1, value = 0.1, step = 0.1)
      )
    ),
    
    mainPanel(
      verbatimTextOutput("results"),
      plotOutput("density_plots")
    )
  )
)




## Define server
**Server Function: server <- function(input, output) {...}

**This defines the server-side logic of the Shiny app, including the simulation process.
Reactive Simulation: run_simulation <- reactive({ ... })

**Defines a reactive expression that runs the simulation whenever input values change.
Initialize Variables:

set.seed(123): Ensures reproducibility of random results.
beta_0_estimates, beta_1_estimates, and beta_2_estimates: Arrays to store the estimated coefficients for each simulation run.
beta_2_estimates is set to NA if procedure2 is chosen, as it doesn’t use beta_2.
Simulation Loop:

**For each iteration (1:input$num_sim), the following steps are executed:
mean_vector and var_covar: Defines mean and covariance matrices for generating correlated random variables x.
x <- mvrnorm(...): Generates correlated variables x1 and x2 based on the correlation (rho) and sample size (n).
epsilon <- rnorm(...): Generates error terms with a specified standard deviation (sigma).
y <- ...: Constructs the response variable using the specified intercept and slopes.
Procedure Logic:

procedure1: Runs a multiple linear regression on both x1 and x2, storing all estimated coefficients.
procedure2: Runs a single linear regression on x1, storing beta_0 and beta_1 estimates only.
procedure3: Runs a multiple regression on both x1 and x2. If the p-value of x2 is below 0.05, it keeps both variables; otherwise, it refits the model using only x1.

server <- function(input, output) {
  
  ### Reactive function to run simulation based on chosen procedure
  run_simulation <- reactive({
    set.seed(123)
    
    # Containers for results
    beta_0_estimates <- numeric(input$num_sim)
    beta_1_estimates <- numeric(input$num_sim)
    beta_2_estimates <- if (input$procedure != "procedure2") numeric(input$num_sim) else NA
    
    for (i in 1:input$num_sim) {
      mean_vector <- c(0, 0)
      var_covar <- matrix(c(1, input$rho, input$rho, 1), nrow = 2)
      
      # Generate random variables for each simulation
      x <- mvrnorm(input$n, mu = mean_vector, Sigma = var_covar)
      epsilon <- rnorm(input$n, mean = 0, sd = input$sigma)
      y <- input$beta_0 + input$beta_1 * x[, 1] + input$beta_2 * x[, 2] + epsilon
      
      # Procedure logic
      if (input$procedure == "procedure1") {
        model <- lm(y ~ x[, 1] + x[, 2])
        beta_0_estimates[i] <- coef(model)[1]
        beta_1_estimates[i] <- coef(model)[2]
        beta_2_estimates[i] <- coef(model)[3]
        
      } else if (input$procedure == "procedure2") {
        model <- lm(y ~ x[, 1])
        beta_0_estimates[i] <- coef(model)[1]
        beta_1_estimates[i] <- coef(model)[2]
        
      } else if (input$procedure == "procedure3") {
        model <- lm(y ~ x[, 1] + x[, 2])
        p_value <- summary(model)$coefficients[3, "Pr(>|t|)"]
        
        # Check that p_value is a single value
        if (length(p_value) == 1 && p_value < 0.05) {
          beta_0_estimates[i] <- coef(model)[1]
          beta_1_estimates[i] <- coef(model)[2]
          beta_2_estimates[i] <- coef(model)[3]
        } else {
          model <- lm(y ~ x[, 1])
          beta_0_estimates[i] <- coef(model)[1]
          beta_1_estimates[i] <- coef(model)[2]
        }
      }
    }
    
    # Return results as a list
    list(
      beta_0_estimates = beta_0_estimates,
      beta_1_estimates = beta_1_estimates,
      beta_2_estimates = if (input$procedure != "procedure2") beta_2_estimates else NA,
      mean_beta_0 = mean(beta_0_estimates),
      mean_beta_1 = mean(beta_1_estimates),
      mean_beta_2 = if (input$procedure != "procedure2") mean(beta_2_estimates, na.rm = TRUE) else NA,
      sd_beta_0 = sd(beta_0_estimates),
      sd_beta_1 = sd(beta_1_estimates),
      sd_beta_2 = if (input$procedure != "procedure2") sd(beta_2_estimates, na.rm = TRUE) else NA,
      last_model = model  # Save the last model for displaying summary
    )
  })
  
  ### Display summary results

**output$results <- renderPrint({ ... }):

**This function generates a print output for display in the UI.
Run the Simulation:

**results <- run_simulation(): Calls the reactive simulation function to get results.
model <- results$last_model: Retrieves the last fitted model from the simulation (if stored).
Print True and Estimated Values:

**cat("True Intercept (beta_0):", input$beta_0, "\n"): Prints the true intercept (user-specified).
cat("Mean of Estimated Intercepts:", results$mean_beta_0, "\n"): Prints the mean of estimated intercepts across all simulations.
cat("Standard Deviation of Estimated Intercepts:", results$sd_beta_0, "\n\n"): Prints the standard deviation of estimated intercepts.
For beta_1 (slope):
Similar cat statements display the mean and standard deviation for beta_1 estimates.
For beta_2 (second slope):
If procedure2 is not selected, displays mean and standard deviation of beta_2 estimates.
Display the Last Model Summary:

**if (!is.null(model)) { print(summary(model)) }: If a model is available, it prints a detailed summary of the last model, including estimates, standard errors, t-values, and p-values for each coefficient.
  
  output$results <- renderPrint({
    results <- run_simulation()
    model <- results$last_model  # Retrieve the last fitted model
    
    cat("True Intercept (beta_0):", input$beta_0, "\n")
    cat("Mean of Estimated Intercepts:", results$mean_beta_0, "\n")
    cat("Standard Deviation of Estimated Intercepts:", results$sd_beta_0, "\n\n")
    
    cat("Mean of Estimated Slope (beta_1):", results$mean_beta_1, "\n")
    cat("Standard Deviation of Estimated Slope (beta_1):", results$sd_beta_1, "\n\n")
    
    if (input$procedure != "procedure2") {
      cat("Mean of Estimated Slope (beta_2):", results$mean_beta_2, "\n")
      cat("Standard Deviation of Estimated Slope (beta_2):", results$sd_beta_2, "\n\n")
    }
    
    # Display the summary of the last model
    
    cat("Display's the summary of the last model from the simulation")
    if (!is.null(model)) {
      print(summary(model))
    }
  })
  
  ### Plot density plots for estimates

**output$density_plots <- renderPlot({ ... }):

**This function creates plots in response to user inputs and displays them in the UI.
Run Simulation:

**results <- run_simulation(): Calls the reactive simulation function to obtain coefficient estimates.
Plot Layout:

**par(mfrow = c(1, ifelse(input$procedure == "procedure2", 2, 3))): Sets up a layout with 2 or 3 columns based on the chosen procedure.
If procedure2 is selected, only two plots (for beta_0 and beta_1) will be shown.
For other procedures, three plots (including beta_2) are shown.
Helper Function for Density Plotting:

**plot_density_safe: A custom function to plot the density of an estimate safely.
Checks if there’s enough variation in the data (more than one unique value) to produce a meaningful density plot.
If yes, plots the density; otherwise, displays a message saying "Insufficient data for density plot."
Density Plots:

**plot_density_safe(results$beta_0_estimates, "blue", ...): Plots the density of intercept (beta_0) estimates in blue.
plot_density_safe(results$beta_1_estimates, "green", ...): Plots the density of beta_1 slope estimates in green.
For beta_2 (if applicable):
If the procedure isn’t procedure2, it plots beta_2 estimates in red.
  
  output$density_plots <- renderPlot({
    results <- run_simulation()
    par(mfrow = c(1, ifelse(input$procedure == "procedure2", 2, 3)))
    
    # Helper function to safely plot densities
    plot_density_safe <- function(data, color, main, xlab) {
      if (length(unique(data[!is.na(data)])) > 1) {
        plot(density(data, na.rm = TRUE), col = color, lwd = 2, main = main, xlab = xlab)
      } else {
        plot.new()
        text(0.5, 0.5, "Insufficient data for density plot", cex = 1.2)
      }
    }
    
    # Plot density for beta_0 and beta_1
    plot_density_safe(results$beta_0_estimates, "blue", "Density of Intercept Estimates", "Estimated beta_0")
    plot_density_safe(results$beta_1_estimates, "green", "Density of Slope Estimates (beta_1)", "Estimated beta_1")
    
    # Plot density for beta_2 if applicable
    if (input$procedure != "procedure2") {
      plot_density_safe(results$beta_2_estimates, "red", "Density of Slope Estimates (beta_2)", "Estimated beta_2")
    }
  })
}

## Run the application 
shinyApp(ui = ui, server = server)
