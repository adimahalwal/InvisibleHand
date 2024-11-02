# InvisibleHand
Monte Carlo simulations 

#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(shiny)
library(tidyverse)
library(MASS)
library(rsconnect)


# Define UI
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

# Define server
server <- function(input, output) {
  
  # Reactive function to run simulation based on chosen procedure
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
  
  # Display summary results
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
  
  # Plot density plots for estimates
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

# Run the application 
shinyApp(ui = ui, server = server)
