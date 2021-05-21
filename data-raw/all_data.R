library(mineCETSA)
#the PI3K file
PI3K1h6h_file <- ms_fileread("210503_1729_PI3K1h6h.txt")
names(PI3K1h6h_file) <- str_replace(names(PI3K1h6h_file), ".x", "1h") #change condition name to fit the treatmentlevel
names(PI3K1h6h_file) <- str_replace(names(PI3K1h6h_file), ".y", "6h")

hitlist_PI3K1h6h <- read.csv("./Data/PI3K/PI3K1h6h_ForHitGeneration_1722_30-03-21_Summary.csv",
                             row.names = 1)
NN_PI3K1h6h <- read.csv("./Data/PI3K/PI3K1h6h_ForHitGeneration_1722_30-03-21_NN.csv",
                        row.names = 1)

#the TNF file
TNF_MOLM1316 <- ms_fileread("210519_1416_TNF_MOLM1316.txt")

hitlist_TNF16 <- read.csv("./Data/TNF/210409_1222_MOLM16_s1_1327_9-04-21_Summary.csv",
                          row.names = 1)
NN_TNF16 <- read.csv("./Data/TNF/210409_1222_MOLM16_s1_1327_9-04-21_NN.csv",
                     row.names = 1)
hitlist_TNF13 <- read.csv("./Data/TNF/210409_1223_MOLM13_s1_1328_9-04-21_Summary.csv",
                          row.names = 1)
NN_TNF13 <- read.csv("./Data/TNF/210409_1223_MOLM13_s1_1328_9-04-21_NN.csv",
                     row.names = 1)
hitlist_TNF <- rbind(hitlist_TNF13, hitlist_TNF16)
NN_TNF <- rbind(NN_TNF13, NN_TNF16)



#the list named drug_data
drug_data <- list("data" = list("PI3K" = PI3K1h6h_file, "TNF" = TNF_MOLM1316),
                  "treat_level" = list("PI3K" = get_treat_level(PI3K1h6h_file), "TNF" = get_treat_level(TNF_MOLM1316)),
                  "hitlist" = list("PI3K" = hitlist_PI3K1h6h, "TNF" = hitlist_TNF),
                  "NN" = list("PI3K" = NN_PI3K1h6h, "TNF" = NN_TNF))

#a graph to print on the tab interactive cell of the app
library(ggplot2)
ev_null_print <-
  ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
  geom_text(x=0.5, y=0.5, label = "Click on a protein \nto see its bar plot !", size = 10) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())


#the interactive cell
library(EBImage)
img2 <- readImage("./data-raw/Animal_Cell.png")

img2@.Data <- aperm(img2@.Data, c(2,1,3))

library(plotly)
fig2 <- plot_ly(type="image", z=img2*255,
                hovertemplate = paste('<br>x: %{x}',
                                      '<br>y: %{y}',
                                      '<extra></extra>')) %>%
  add_trace(type = "scatter",
            mode = "markers",
            x = loca_orga$x[c(5, 240)], y = loca_orga$y[c(5, 240)],
            customdata = c("P85037", "Q9Y3U8"),
            text = c("P85037", "Q9Y3U8"),
            hovertext = loca_orga$organelle[c(5, 240)],
            marker =
              list(color = c("red", "blue"),
                   size = 10),
            hovertemplate = paste('Protein: %{text}',
                                  '<br>Organelle: %{hovertext}',
                                  '<extra></extra>'),
            showlegend = FALSE
  ) %>%
  layout(xaxis = list(visible = FALSE),
         yaxis = list(visible = FALSE))


usethis::use_data(PI3K1h6h_file, overwrite = TRUE)
usethis::use_data(hitlist_PI3K1h6h, overwrite = TRUE)
usethis::use_data(NN_PI3K1h6h, overwrite = TRUE)
usethis::use_data(TNF_MOLM1316, overwrite = TRUE)
usethis::use_data(hitlist_TNF, overwrite = TRUE)
usethis::use_data(NN_TNF, overwrite = TRUE)
usethis::use_data(drug_data, overwrite = TRUE)
usethis::use_data(ev_null_print, overwrite = TRUE)
usethis::use_data(img2, overwrite = TRUE)
usethis::use_data(fig2, overwrite = TRUE)


