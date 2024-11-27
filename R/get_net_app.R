#' get_net_app
#'
#' Function to get the network from specific proteins. This function is totally based on the function plot_network from
#' the STRINGdb package, but it now return the plot (using ggdraw) or return an interactive image (using plotly).
#'
#' @param string_ids The STRING ids of the proteins you want to plot the network
#' @param payload_id An identifier of payload data on the STRING server (see method post_payload
#'                   for additional information)
#' @param required_score A threshold on the score that overrides the default score_threshold, that we use
#'                       only for the picture
#' @param inter A logical to tell if you want to return an interactive image with plotly
#' @param add_link Whether you want to generate and add a short link to the
#'                 relative page in STRING. As default this option is active but we suggest to deactivate it in case
#'                 one is generating many images (e.g. in a loop). Deactivating
#'                 this option avoids to generate and store a lot of short-urls on our server.
#' @param network_flavor Specify the flavor of the network
#'                       ("evidence", "confidence" or "actions". default "evidence").
#' @param add_summary Parameter to specify whether you want to add a summary text to the picture.
#'                    This summary includes a p-value and the number of proteins/interactions.
#'
#' @return The network
#'
#' @export
#'
#' @seealso \code{\link{STRINGdb}}
#'

get_net_app <- function(string_ids, payload_id = NULL,
                        required_score = NULL, inter = FALSE,
                        add_link = TRUE, network_flavor = "evidence",
                        add_summary = TRUE){
  if (is.null(required_score))
    required_score = string_db$score_threshold
  img_t <- string_db$get_png(string_ids, payload_id = payload_id,
                              required_score = required_score,
                              network_flavor = network_flavor)

  if(!inter){
    if (!is.null(img_t)) {
      p <- ggdraw() + draw_image(img_t, width = 1.5, x = -0.25)

      if (add_summary)
        p <- p + draw_label(
          string_db$get_summary(string_ids, required_score),
          x = 0.25, y = 0.05, size = 9,
        )
    }
  }
  else{
    if (!is.null(img_t)) {
      p <- plot_ly(type="image", z = img_t*255)


      if (add_summary)
        p <- p %>%
          layout(title = list(text = string_db$get_summary(string_ids, required_score),
                              font = list(size = 11,
                                          color = "black")
          ),
          margin = list(t = 90, r = 5, l = 5, b = 5)
          )

      p <- p %>%
        layout(xaxis = list(visible = FALSE),
               yaxis = list(visible = FALSE))
    }
  }
  return(p)
}
