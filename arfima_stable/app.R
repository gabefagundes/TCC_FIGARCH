#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(stabledist)
library(Cairo)
library(plotly)
library(ggplot2)
options(shiny.usecairo = T)

source("ARFIMA 2.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Teste Variáveis Alpha Estáveis ARFIMA"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput('n',
                         "Número de observações",
                         value = 500),
            numericInput('iter',
                         "Número de Interações",
                         value = 50),
            sliderInput("d",
                        "d (deve ser menor que 1 - 1/alpha)",
                        min = -0.5,
                        max = 0.5,
                        value = 0.3,
                        step = 0.05),
            sliderInput("alpha",
                        "Alpha",
                        min = 0,
                        max = 2,
                        value = 2,
                        step = 0.1),
            sliderInput("beta",
                        "Beta",
                        min = -1,
                        max = 1,
                        step = 0.1,
                        value = 0),
            numericInput('gamma',
                         "Gamma",
                         value = 1), 
            numericInput('delta',
                         "Delta",
                         value = 0),
            width = 2
        ),
        

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("distPlot",
                      height = '800px'),
           width = 10
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlotly({
        
        y = arfima_proc(0, input$d, 0, 
                        n = input$n, innov = rstable, iter = input$iter,
                        alpha = input$alpha, beta = input$beta, 
                        gamma = input$gamma, delta = input$delta)
        
        ggplotly(
            ggplot()+
                geom_line(aes(x = 1:input$n, y = y)) +
                theme_minimal() +
                xlab('t')+
                ylab('y_t')+
                ggtitle("Teste para processo FI(d) com inovações alpha-estáveis") +
                scale_x_continuous(breaks = seq(0, input$n, length = 21))
        )
        
    })
    
    output$max_d = renderText({
        1 - 1/input$alpha
    }
    )
        
}

# Run the application 
shinyApp(ui = ui, server = server)
