#PP-plot function
.ppPlot <- function(data,ylim=NULL){
pi0 <- sum(data$decoy)/sum(!data$decoy)
ppPlot <- ggplot()  +
  geom_abline(slope = pi0,color = 'black') +
  labs(title = 'PP plot') +
  coord_cartesian(xlim = c(0,1), ylim = ylim, expand = TRUE) +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.2)),
    axis.text = element_text(size = rel(1.2)),
    axis.title.y = element_text(angle = 0))
x <- data$score[!data$decoy]
Ft <- ecdf(x)
Fd <- ecdf(data$score[data$decoy])
df <- data_frame(Fdp = Fd(x), Ftp = Ft(x))
ylimHlp <- mean(Fd(x)==1)

return(list(ppPlot=ppPlot + geom_point(data = df,aes(Fdp,Ftp),color = 'dark gray'),
yLim=c(0,(1-ylimHlp)),
Fd=Fd,
Ft=Ft))
}

#Function to open pop-up window where variables can be chosen
.select <- function(object,decoy=NULL,score=NULL,log=TRUE,nBins=50) {
    require(shiny)
    require(ggplot2)
    out <- list(selDecoy=decoy,selScore=score,log=log,nBins=nBins)
    fv <- colnames(object)
    on.exit(return(out))
    ui <- shiny::fluidPage(
      title = 'Evaluate Decoys',
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::actionButton("stop", "Stop app"),
          checkboxInput("log", "-log10 transform variable?", out$log),
          shiny::selectInput("decoyVar", "Select Decoy",
                                    as.list(fv), selected = out$selDecoy),
          shiny::selectInput("scoreVar", "Select Score",
                                    as.list(fv), selected = out$selScore),
          shiny::numericInput("nBins", "Number of bins in histogram", value = out$nBins, min = 2, max = 1000)
        ),
        shiny::mainPanel(
        tabsetPanel(type = "tabs",
                 tabPanel("Info", HTML(
                  paste(
                  h4("Use the left panel to"),
                  '<ul>
                    <li>Indicate if the scores should be log-transformed;</li>
                    <li>Select a <i>Decoy</i> variable</li>
                    <li>Select a <i>Score</i> variable</li>
                    <li>Select the number of <i>Bins</i> for the histogram </li>
                  </ul>',
                  h4("The following conditions must be met"),
                  '<ul>
                  <li> Decoy variable is a boolean that indicates if the score belongs to a target or a decoy</li>
                  <li> Score variable contains the scores of the search engine, which have to be continuous (larger scores are assumed to be better. E-values are typically -log10(e-value) transformed.)</li>
                  <li> You can verify the selected data and settings using the table below, and, the plots in the histogram and PP-Plot tabs.</li>
                  </ul>',
                  "Press button <b>Stop app</b> to return to R, and, to generate and store the plot in your R session."
                  )
                  ),
                  h4("Data table with selected variables"),
                  shiny::dataTableOutput('fd')
                  ),
                  tabPanel("Histogram",shiny::plotOutput('histPlot',
                                          dblclick = "hist_dblclick",
                                          brush = brushOpts(id = "hist_brush",resetOnNew = TRUE)),
                                          p("Brush and double-click on a selected area to zoom in on x-coordinate.
                                          Double click outside a selected area to zoom out")
                                          ),
                  tabPanel("PP-Plot",shiny::plotOutput('PPplot',
                                          dblclick = "PPplot_dblclick",
                                          brush = brushOpts(id = "PPplot_brush",resetOnNew = TRUE)),
                                          p("Brush and double-click on a selected area to zoom in on y-coordinate.
                                          Double click outside a selected area to zoom out")
                                          )
                  )

        )
        )
    )

    server <- function(input, output) {
      # Terminate app
      shiny::observeEvent(input$stop, {
        shiny::stopApp(returnValue = out)
      })

      # Make Data Table
      output$fd <- shiny::renderDataTable({
        out$selDecoy <<- input$decoyVar
        out$selScore <<- input$scoreVar
        out$log <<- input$log
        out$nBins <<- input$nBins
        object[, c(out$selDecoy,out$selScore), drop = FALSE]
      })

      # Select Data
      data <- reactive({
              data <- object[,c(input$decoyVar,input$scoreVar)]
              data <- na.exclude(data)
              names(data)<-c("decoy","score")
              data$score<-as.double(data$score)
              if (input$log) data$score<--log10(data$score)
              data
      })


      ############
      # Histogram#
      ############

          output$histPlot <- renderPlot({
          binwidth <- diff(range(data()$score,na.rm=TRUE))/input$nBins
          basePlot <- ggplot(data(),
             aes(score, fill = decoy, col=I("black")))
          if(is.logical(data()$decoy)&is.numeric(data()$score))
          basePlot +
             geom_histogram(alpha = 0.5, bins=input$nBins, position = 'identity') +
             labs(x = 'Score', y = 'Counts' ,title = paste0('Histogram of targets and decoys\n')) +
             coord_cartesian(xlim = histRanges$x, ylim = histRanges$y, expand = TRUE) +
             theme_bw() +
             theme(
               plot.title = element_text(size = rel(1.5)),
               axis.title = element_text(size = rel(1.2)),
               axis.text = element_text(size = rel(1.2)),
               axis.title.y = element_text(angle = 0))
              })

              histRanges <- reactiveValues(x = NULL,y=NULL)

              observeEvent(input$hist_dblclick, {
                  brush <- input$hist_brush
                  if (!is.null(brush))
                    histRanges$x <- c(brush$xmin, brush$xmax) else
                    histRanges$x <- NULL
                })

      #########

      #########
      #PP-plot#
      #########
            output$PPplot <- renderPlot({
            basePlot<-ggplot()
            if(is.logical(data()$decoy)&is.numeric(data()$score))
            basePlot<-.ppPlot(data(),PPplotRanges$y)[[1]]
            basePlot
            })
            PPplotRanges <- reactiveValues(x = c(0,1),y=c(0,1))

            observeEvent(input$PPplot_dblclick, {
                brush <- input$PPplot_brush
                if (!is.null(brush))
                  PPplotRanges$y <- c(0, brush$ymax) else
                  PPplotRanges$y <- c(0,1)
              })
      ##########
    }
    app <- list(ui=ui, server=server)
    shiny::runApp(app)
}


##' Function to evaluate target/decoy searches
##'
##' @description Evaluation of target/decoy approach
##'
##' @param object An \code{mzIDfile}, \code{mzRfile} or \code{dataframe}.
##'
##' @param decoy A \code{boolean} (default is \code{NULL}) to indicate whether one has a target or a decoy.
##'
##' @param score A \code{numeric} variable (default is \code{NULL}), containing scores of the engine.
##'
##' @param log A \code{logical} (default is \code{TRUE}) indicating whether the score
##' variable should be -log10-transformed.
##'
##' @return PP-plots and histograms to check if the decoys are a good simulations of the bad target hits.
##'
##' @example
##'
##' library(mzID)
##' exampleFiles <- list.files(system.file('extdata', package = 'mzID'),
##'                           pattern = '*.mzid', full.names = TRUE)
##' mzIDexample<-mzID(exampleFiles[[2]])
##' decoyPlots<-evalTargetDecoys(object=mzIDexample,decoy="isdecoy",score="x\\!tandem:expect",log10=TRUE,nBins=30)
##' decoyPlots$together
##'
##' #Individual plots can also be plotted
##' names(decoyPlots)
##' decoyPlots$histogramZoom
##'
##' # If you do not know the name of the score and/or decoy variable,
##' # or if you want to evaluate how many bins you want to use in the histogram,
##' # or if -log10 transformation is needed you can launch a shiny app by only
##' # specifying the mzID object
##'
##' # evalTargetDecoys(mzIDexample)
##' @author
##' @export

evalTargetDecoys <- function(object,
                             decoy=NULL,
                             score=NULL,
                             log10=NULL,
                             nBins=50)
  {
  # require some packages
  require(mzID)
  require(mzR)
  require(dplyr)
  require(ggplot2)
  require(graphics)
  require(ggpubr)
  require(shiny)

  # check object class
  if (class(object) == "mzID") {
    df <- flatten(object)
  } else if(class(object) == "data.frame") {
    df <- object
  } else if(class(object) == "mzRident") {
    df <- cbind(psms(object),score(object)[,-1])
  } else {"object should be of the class mzID, mzRident or dataframe"}

  # if one or more arguments in the function are missing, .select() is called, a pop-up window opens and the variables must be chosen manually.
  if(missing(decoy)|missing(score)|missing(log10)) {
  if (missing(log10)) log10 <-TRUE
  out <- .select(df,decoy,score,log10,nBins)
  decoy <- out$selDecoy #categorical
  score <- out$selScore #continu
  log10 <- out$log
  nBins <- out$nBins
  }

  # subset of dataframe
  data <- df[, c(decoy, score)]
  names(data) <- c("decoy", "score")
  data <- na.exclude(data)
  data$score<-as.double(data$score)

  # if variable 'score' is a character, change to continuous
  if(class(data$score) == "character"){
    data$score <- as.numeric(as.character(data$score))
  }

  # check whether the selected variables are of the correct class
  if(class(data$score)!="numeric") stop('Score is not numeric')
  if(class(data$decoy)!="logical") stop('Decoy is not logical')

  # perform log10-transformation on variable 'score' if so indicated
  if(log10) {data$score <- -log10(as.numeric(data$score))}

#############
### Plots ###
#############

  # create PP-plot
  p1<-.ppPlot(data)

  # create histogram
  p2 <- ggplot(data, aes(score, fill = decoy, col=I("black"))) +
  geom_histogram(alpha = 0.5, bins=nBins, position = 'identity') +
    labs(x = head(score), y = "",
         title = 'Histogram of targets and decoys') +
    theme_bw() +
    theme(
      plot.title = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.2)),
      axis.text = element_text(size = rel(1.2)),
      axis.title.y = element_text(angle = 0))


  out<-list(ppPlot=p1$ppPlot,
      histogram=p2,
      ppPlotZoom=p1$ppPlot+ylim(p1$yLim[1],p1$yLim[2]),
      histogramZoom=p2+coord_cartesian(xlim = c(min(data$score),max(data$score[data$decoy])),expand=TRUE)
  )
  out$together<-ggarrange(plotlist=out,ncol=2,nrow=2)
  return(out)
}
