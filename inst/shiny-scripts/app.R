# This example is adapted from
# Grolemund, G. (2015). Learn Shiny - Video Tutorials. URL:https://shiny.rstudio.com/tutorial/

library(shinyalert)
library(shiny)

# Define UI
ui <- fluidPage(

  # Change title
  titlePanel(tags$b("orthoVisualizer:"),"Quantify, visualize, and annotate orthologs in DNAseq"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("There are 4 functions in this package. 1) annotateSeq: Prints DNAseq
             in lines of 50 nucleotides and underlines corresponding orthologs
             found. Returns invisible null. 2) quantSeq: Returns tibble of
             5 columns (id, header, number of ortholog motifs found,
             length of sequence, frequency ratio of motif) for each sequence in
             the fasta file of interest. 3) freqSeq: Displays and returns
             occurrence of ortholog gene/motif subsequence for each DNAseq in
             target fasta file. 4) freqRatioSeq: Displays and returns frequency
             ratio of ortholog gene/motif subsequence for each DNAseq in target
             fasta file. See vignette under freqRatioSeq for details on
             frequency ratio."),
      # br() element to introduce extra vertical spacing ----
      br(),

      tags$b("Description: This Shiny App is part of the orthoVisualizer. It
             permits to identify occurrences of ortholog gene/motif subsequence
             for each DNAseq in the target fasta file to annotate, quantify, and
             visualize abundance of such subsequences. For more details, see
             vignette."),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),

      # input
      tags$p("Instructions: Below, enter or select values required to perform the analysis.
             Default values are shown. Then press 'Run'. Navigate through
             the different tabs to the right to explore the results."),

      # br() element to introduce extra vertical spacing ----
      br(),

      # input
      shinyalert::useShinyalert(force = TRUE),  # Set up shinyalert

      textInput(inputId = "fastaPath",
                label = "Enter a character string representing the filepath for the
                target DNAseq fasta file.",
                system.file("extdata", "DNA_seq.fasta",
                            package = "orthoVisualizer")),
      textInput(inputId = "pattern",
                label = "Enter a character string representing the ortholog
                gene/motif subsequence to search for", "CGCG"),

      # br() element to introduce extra vertical spacing ----
      br(),

      # actionButton
      actionButton(inputId = "button2",
                   label = "Run"),

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Annotation of Orthologs in DNASeq",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h4("Ortholog Annotation:"),
                           br(),
                           verbatimTextOutput("annotate")),
                  tabPanel("Quantification of Ortholog Gene/Motif in DNASeqs",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Quantification:"),
                           br(),
                           verbatimTextOutput("quantify")),
                  tabPanel("Plot of Ortholog Frequency (Number of Occurrence)",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Plot Ortholog Frequency (Number of Occurrence):"),
                           br(),
                           br(),
                           plotOutput("DNAfreqPlot")),
                  tabPanel("Plot of Ortholog Frequency Ratio",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Plot Ortholog Frequency Ratio:"),
                           br(),
                           br(),
                           plotOutput("DNAfreqRatioPlot"))
      )

    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {

  # Quantify
  startquantification <- eventReactive(eventExpr = input$button2, {
    orthoVisualizer::quantSeq(
      fastaPath = input$fastaPath,
      pattern = input$pattern)
  })


  # Text output for quantification of Ortholog appearance in DNASeqs
  output$quantify <- renderPrint({
    if (! is.null(startquantification))
      startquantification()
  })

  # Annotation of orthologs
  output$annotate <- renderPrint({
    if (! is.null(startquantification))
      orthoVisualizer::annotateSeq(fastaPath =  input$fastaPath,
                                   pattern = input$pattern)
  })


  # Plotting frequency of pattern in DNAseqs in fastafile
  output$DNAfreqPlot <- renderPlot({
    if (! is.null(startquantification)) {
      orthoVisualizer::freqSeq(fastaPath = input$fastaPath,
                               pattern = input$pattern)
    }
  })

  # Plotting frequency of pattern in DNAseqs in fastafile
  output$DNAfreqRatioPlot <- renderPlot({
    if (! is.null(startquantification)) {
      orthoVisualizer::freqRatioSeq(fastaPath = input$fastaPath,
                                    pattern = input$pattern)
    }
  })


}

# Create Shiny app ----
shiny::shinyApp(ui, server)

# [END]
