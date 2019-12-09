library(shiny)
library(ggplot2)
library(VennDiagram)

se <- function(x) sqrt(var(x)/length(x))

ui <- fluidPage(

    # Application title
    titlePanel("Biostatisztika GR"),

    navlistPanel(
        "#1 Becslés",
        tabPanel(
            title="1. Pontbecslés",
            sidebarLayout(
                sidebarPanel(
                    radioButtons('sampleSize', label=h3('Mintanagyság'),
                                 choices=list(10,100,1000),
                                 selected=10),
                ),
                mainPanel(
                    plotOutput("popPlot")
                )
            )
        ),
        tabPanel(
            title="2. Intervallumbecslés",
            mainPanel(
                plotOutput("intervalplot")
            )
        ),
        "#2 Valószínűségszámítás",
        tabPanel(
            title="1. Valószínűségszámítás",
                mainPanel(
                    plotOutput("vennplot")
            )
        ),
        tabPanel(
            title="2. Valószínűségi változók",
            sidebarLayout(
                sidebarPanel(
                    radioButtons('dist', label=h3('Eloszlás'),
                                 choices=c("Egyenletes" = "uniform",
                                           "Binomiális" = "binom",
                                           "Poission" = "pois",
                                           "Normál" = "normal")),
                    h3("Egyenletes eloszlás paraméterei"),
                    sliderInput("min", label="min", value=0, min=0,max=10),
                    sliderInput("max", label="max", value=1, min=1,max=20),
                    h3("Binomiális eloszlás paraméterei"),
                    sliderInput("size", label="n", value=10, min=1,max=50),
                    sliderInput("prob", label="p", value=0.5, min=0.1,max=0.9, step=0.1),
                    h3("Poisson eloszlás paraméterei"),
                    sliderInput("lambda", label="λ", value=4, min=1,max=50),
                    h3("Normál eloszlás paraméterei"),
                    sliderInput("mean", label="μ", value=0, min=-10,max=10),
                    sliderInput("sd", label="σ", value=1, min=1,max=10),
                ),
                mainPanel(
                    plotOutput("distplot")
                )
            )
        ),
        tabPanel(
            title="3. Központi határeloszlás tétele",
            mainPanel(uiOutput("clt"))
        ),
        "#3 Epidemiológia",
        tabPanel(
            title="1. Epidemiológia",
            mainPanel(
                imageOutput("epidemiology")
            )
        ),
        tabPanel(
            title="2. Teljes valószínűség tétele",
            mainPanel(
                imageOutput("tvt"),
                withMathJax(),
                helpText("A teljes valószínűség tétele:
                         $$P(A) = \\sum P(A \\cap B_x) = \\sum P(A|B_x) \\ast P(B_x)$$")
            )
        ),
        "#4 Hipotézisvizsgálat",
        tabPanel(
            title="1. Hipotézisvizsgálat",
            mainPanel(
                imageOutput("hipottable")
            )
        )
    )
)

# Define server logic
server <- function(input, output, session) {
    
    output$popPlot <- renderPlot({
        x=seq(-4, 4, length=10000)
        y=dnorm(x)
        sample_x=sample(x, size=input$sampleSize, replace=FALSE)
        sample_se=se(sample_x)
        sample_mean=mean(sample_x)
        plot(x, y, type="l")
        abline(v=0)
        points(sample_x, sample_x*0)
        abline(v=sample_mean, col="red")
        abline(v=sample_mean - sample_se, col="blue")
        abline(v=sample_mean + sample_se, col="blue")
        mtext("")
        mtext(paste("Mintaátlag: ", format(sample_mean, digits=1), "±", format(sample_se, digits=2)))
    })
    
    output$intervalplot <- renderPlot({
        samsize <- 100 
        replicates <- 50
        pval <- .05
        
        samples <- replicate(replicates, rnorm(samsize))
        confint <- t(apply(samples, 2, function(x)
            c(mean(x)-qt(1-pval/2, 
                         df=samsize-1)*sd(x)/sqrt(samsize), 
              mean(x)+qt(1-pval/2, df=samsize-1)*sd(x)/sqrt(samsize))))
        
        # Simple plot
        plot(c(0, 0), c(1, replicates), col="black", typ="l",
             ylab="Samples",
             xlab="Confidence Interval")
        segments(confint[,1], 1:replicates, confint[,2], 1:replicates)
        
        # Use red if mean outside interval
        outside <- ifelse(confint[,1]>0 | confint[,2]<0, 2, 1)
        plot(c(0, 0), c(1, replicates), col="black", typ="l",
             ylab="Minták",
             xlab="Konfidencia-intervallumok")
        segments(confint[,1], 1:replicates, confint[,2], 1:replicates,
                 col=outside)
    })
    
    output$vennplot <- renderPlot({
        draw.pairwise.venn(area1 = 0.5, area2 = 0.3, cross.area = 0.1, cat.pos=0, category=c(
            "Mozgásszervi probléma", "Érzékszervi probléma"))
    })
    
    output$distplot <- renderPlot({
        dist <- switch(input$dist,
                       uniform= runif(500, min=input$min, max=input$max),
                       binom=rbinom(500, size=input$size, prob=input$prob),
                       pois=rpois(500, lambda=input$lambda),
                       normal=rnorm(500, mean=input$mean, sd=input$sd),
                       runif)
        hist(dist)
    })
    
    output$clt <- renderUI({
        tagList("URL: ", a("Shiny app a központi határeloszlás tételének vizualizálására", 
                           href="https://gallery.shinyapps.io/CLT_mean/"))
    })
    
    output$epidemiology <- renderImage({
        return(list(
            src = "images/epidem.png",
            filetype = "image/png",
            alt = "Az epidemiológus fürdőkádja",
            height = "600px"
        ))
    }, deleteFile = FALSE)
    
    output$tvt <- renderImage({
        return(list(
            src = "images/tvt.png",
            filetype = "image/png",
            alt = "A teljes valószínűség tételének vizualizációja",
            height="400px"
            
        ))
    }, deleteFile = FALSE)
    
    output$hipottable <- renderImage({
        return(list(
            src = "images/hypot.png",
            filetype = "image/png",
            alt = "Hipotézisvizsgálat igazságtábla",
            width = "900px"
        ))
    }, deleteFile = FALSE)
    
    session$onSessionEnded(function() {
        stopApp()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
