#' hit_plotcell
#'
#' Function to return the 'cell plot' of the proteins according to their category.
#'
#' @param data The output from hit_for_cell function
#' @param tit The title of the plot
#' @param cond The treatments you want to keep. If NULL, will take it all
#' @param cat_col_list A list containing the color for each category (CC, CN, NC, ND, NN)
#'
#' @return An interactive plot
#'
#' @export
#'
#' @seealso \code{\link{hit_for_cell}}


hit_plotcell <- function(data, tit = "PI3K data in the cell",
                         cond = NULL, cat_col_list = list("CC" = "red", "CN" = "lightblue",
                                                          "NC" = "yellow", "ND" = "#747474",
                                                          "NN" = "#CCCCCC")){
  df <- data
  if(!is.null(cond)){
    df <- df[which(!is.na(match(df$treatment, cond))),]
  }

  df <- df[order(df$nb_location),]
  cat_color <- levels(df$category)
  cat_color <- as.character(cat_col_list[cat_color])

  p <- plot_ly(type="image", z=img2*255,
               colors = cat_color,
               symbols = c("circle", "x", "triangle-down", "square"),
               hoverinfo = "text", source = "M") %>%
    add_trace(type = "scatter", mode = "markers",
              data = df,
              customdata = ~id,
              x = ~x, y = ~y,
              color = ~category, symbol = ~nb_location,
              text = ~txt,
              marker = list(size = 8)
    ) %>%
    layout(xaxis = list(visible = FALSE),
           yaxis = list(visible = FALSE),
           title = list(text = tit,
                        font = list(size = 25,
                                    family = "Times New Roman")
           ),
           margin = list(t = 90, b = 90))

  return(p)
}
