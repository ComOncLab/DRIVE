library(shiny)
library(shinyWidgets)
library(bslib)

ui <- page_sidebar(
  title = "现代美观的文本输出",
  theme = bs_theme(version = 5, bootswatch = "minty"),
  
  sidebar = sidebar(
    numericInput("value", "输入数值:", 100, min = 1, max = 1000),
    textAreaInput("long_text", "输入长文本:", 
                  "这是一段示例文本，用于展示多行文本的输出效果。", 
                  rows = 3),
    prettySwitch("show_details", "显示详细信息", TRUE)
  ),
  
  card(
    card_header("漂亮的卡片展示"),
    htmlOutput("styled_output"),
    value_box(
      title = "统计结果",
      value = textOutput("summary_value"),
      showcase = icon("chart-bar")
    )
  ),
  
  card(
    card_header("进度和状态"),
    progressBar("pb", value = 50, status = "info", striped = TRUE),
    br(),
    alert(
      status = "info",
      icon = icon("info-circle"),
      "这里是提示信息"
    )
  ),
  
  card(
    card_header("可折叠面板"),
    accordion(
      accordion_panel(
        "详细信息",
        icon = icon("magnifying-glass"),
        htmlOutput("details_output")
      )
    )
  )
)

server <- function(input, output) {
  
  output$styled_output <- renderText({
    paste0(
      '<div class="alert alert-primary" role="alert">',
      '<h5><i class="bi bi-info-circle"></i> 分析结果</h5>',
      '<hr>',
      '<div class="row">',
      '<div class="col-md-6">',
      '<p><strong>输入数值:</strong> <span class="badge bg-primary">', 
      input$value, '</span></p>',
      '</div>',
      '<div class="col-md-6">',
      '<p><strong>平方值:</strong> <span class="badge bg-success">', 
      input$value^2, '</span></p>',
      '</div>',
      '</div>',
      '</div>'
    )
  })
  
  output$summary_value <- renderText({
    paste("平均值:", mean(c(input$value, input$value^2)))
  })
  
  output$details_output <- renderText({
    paste0(
      '<div class="p-3 bg-light rounded">',
      '<h6>详细分析:</h6>',
      '<ul>',
      '<li>输入值: ', input$value, '</li>',
      '<li>平方值: ', input$value^2, '</li>',
      '<li>平方根: ', round(sqrt(input$value), 2), '</li>',
      '<li>对数: ', round(log(input$value), 2), '</li>',
      '</ul>',
      '<p class="mb-0"><small><i>最后更新: ', Sys.time(), '</i></small></p>',
      '</div>'
    )
  })
  
  observe({
    updateProgressBar(
      session = getDefaultReactiveDomain(),
      id = "pb",
      value = input$value %% 100
    )
  })
}

shinyApp(ui, server)