# Install required packages if not already installed
# install.packages(c("shiny", "readxl", "ggplot2", "factoextra", "dplyr"))

library(shiny)
library(readxl)
library(ggplot2)
library(factoextra)
library(dplyr)

# Define UI for app
ui <- fluidPage(
  titlePanel("PCA App"),
  
  # Sidebar layout for input and output
  sidebarLayout(
    sidebarPanel(
      # Input: Upload Excel file
      fileInput("file", "Upload Excel File", accept = c(".xlsx")),
      
      # Input: Select numeric columns for PCA
      uiOutput("numeric_cols"),
      
      # Input: Select PC axes to plot
      selectInput("pc_x", "Select PC for X-axis", choices = NULL),
      selectInput("pc_y", "Select PC for Y-axis", choices = NULL),
      
      # Input: Select character/factor column for annotation
      uiOutput("factor_col")
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      # Output: PCA plot
      plotOutput("pca_plot")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive function to load and read the Excel file
  data <- reactive({
    req(input$file)
    read_excel(input$file$datapath)
  })
  
  # Reactive function to segregate columns into numeric and character
  observe({
    req(data())
    df <- data()
    
    # Detect numeric and character columns
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    factor_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    
    # Update UI for numeric columns
    output$numeric_cols <- renderUI({
      checkboxGroupInput("selected_numeric_cols", "Select Numeric Columns for PCA",
                         choices = numeric_cols, selected = numeric_cols)
    })
    
    # Update UI for character/factor column selection
    output$factor_col <- renderUI({
      selectInput("selected_factor_col", "Select Character/Factor Column for Annotation",
                  choices = c("None", factor_cols), selected = "None")
    })
  })
  
  # Reactive function to perform PCA on selected numeric columns
  pca_result <- reactive({
    req(input$selected_numeric_cols)
    df <- data()
    
    # Subset selected numeric columns
    selected_data <- df %>% select(all_of(input$selected_numeric_cols))
    
    # Perform PCA
    prcomp(selected_data, scale. = TRUE)
  })
  
  # Update choices for PC axes once PCA is performed
  observe({
    req(pca_result())
    pcs <- colnames(pca_result()$x)
    
    updateSelectInput(session, "pc_x", choices = pcs, selected = pcs[1])
    updateSelectInput(session, "pc_y", choices = pcs, selected = pcs[2])
  })
  
  # Plot PCA plot with selected PCs and annotation
  output$pca_plot <- renderPlot({
    req(pca_result(), input$pc_x, input$pc_y)
    
    # Extract PCA results and create a data frame
    pca_df <- as.data.frame(pca_result()$x)
    df <- data()
    
    # Check if the user selected a character/factor column for annotation
    if (input$selected_factor_col != "None") {
      pca_df$annotate <- df[[input$selected_factor_col]]
    } else {
      pca_df$annotate <- NULL
    }
    
    # Plot PCA with ggplot
    p <- ggplot(pca_df, aes_string(x = input$pc_x, y = input$pc_y)) +
      geom_point(aes(color = annotate), size = 3) +
      labs(title = "PCA Plot", x = input$pc_x, y = input$pc_y) +
      theme_minimal()
    
    # Return the plot
    print(p)
  })
}

# Run the app
shinyApp(ui = ui, server = server)
