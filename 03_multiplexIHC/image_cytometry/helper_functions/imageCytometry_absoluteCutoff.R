#' @title imageCytometry_absoluteCutoff
#'
#' @description Perform automated image cytometric analysis with reference to master table. It can plot and sort cells using following gates :\iteminze{
#' \item qunadrantGate gate based on X and Y cutoff
#' \item quantileGate gate based on single marker cutoff
#' \item rectangularGate gate based on X and Y rectangular cufoff
#' \item tailGate gate to obtain tail gate (+/-)
#' \item hlineGate provide X cutoff
#' \item vlineGate provide Y cutoff
#' \item boundaryGate gate to cover all cells; no need to provide cutoff
#' \item singletGate gate based on identification of singlet cells
#' \item autoGate gate based on automatic clusters
#' \item # no gate Select complete population
#' }
#'
#' @param data Image cytometry data from the CellProfiler in csv format.
#' @param master.table characteristics master table
#' @param K Number of clusters to use fro autoGating. (dafault is 1)
#' @param eps refer dbscan() function in dbscan package.(dafault is 0.5)
#'
#' @imoort tidyverse
#' @importFrom dichromat colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ks kde
#' @importFrom sp point.in.polygon
#' @importFrom MASS rlm
#' @importFrom scales percent
#' @importFrom plyr ldply
#' @impoerFrom tools toTitleCase
#' @importFrom dbscan dbscan
#' @importFrom  patchwork wrap_plots
#' @import utils
#' @import stats
#' @import graphics
#' @import grDevices
#'
#' @return ImageCytometry provied following results.
#'   The object contains the following:\itemize{
#'     \item plot facs plots of the sample
#'     \item cordinate.data cordiante data with phenotype map
#'     \item summary.results summary of distribution
#'   }
#'
#' @examples
#' imageCytometry(data=data, master.table 0 master.table)
#' @export
imageCytometry_absoluteCutoff <- function(data = data,
                           master.table = master.table,
                           K = 1,
                           eps = 0.5) {

  ## function check
  stopifnot(inherits(data, "data.frame"))
  validObject(data)

  stopifnot(inherits(master.table, "data.frame"))
  validObject(master.table)
  

## master table
  master.table <- master.table %>%
    mutate_all(na_if, "")
  ## create datalist per applied gate
  data.list <- list()
  facs.plot <- list()

  data[["phenotypes"]] <- NA

  data.list[[toupper("Root")]] <- data

  for (mt in 1:nrow(master.table)) {
    ## color palette
    colfunc <- colorRampPalette(rev(brewer.pal(10, 'Spectral')))

    ## ----------------------------------------------------------
    ## root gate
    ## ----------------------------------------------------------

    if (master.table[mt, "appliedGate"] %in% "Root") {
      ## retrieve data for applied gate
      fs <- data.list[[toupper(master.table[mt, "appliedGate"])]]

      if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

        ### workaround for error in counter plots
        drop.cols <- c("IMAGENUMBER", "phenotypes")

        counter.errors <- fs %>%
          .[, setdiff(names(.), drop.cols)] %>%
          gather("variable", "value") %>%
          group_by(variable) %>%
          summarise(
            q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            q75 = quantile(value, probs = 0.75, na.rm = TRUE)
          ) %>%
          mutate(diff = q75 - q25) %>%
          filter(diff == 0)

        if (dim(counter.errors)[1] == 0) {
          ## select x,y
          xy <- toupper(master.table[mt, "markers"])
          xy <- gsub(" ", "", unlist(strsplit(xy, ",")))
          x <- xy[1]
          y <- xy[2]

          ## ----------------------------------------------------------
          ## Quandrant gate
          ## This gate is based on quandrant cutoff for X and Y (need
          ## to provide single cutoff for X and Y)
          ## ----------------------------------------------------------

          if (master.table[mt, "gateType"] %in% "quadrantGate") {
            ## gate cutoff
            gateCutoff <- master.table[mt, "cutoff"]
            gateCutoff <-
              as.numeric(gsub(" ", "", unlist(strsplit(
                gateCutoff, ","
              ))))

            ##get quantile cutoff
            x_quant <- gateCutoff[1]
            y_quant <- gateCutoff[2]

            ## gate to choose
            sign <-
              gsub(" ", "", master.table[mt, "markerDistribution"])

            if (nchar(sign) != 3) {
              print(
                paste0(
                  "please provide binary marker distriubtion for ",
                  master.table[mt, "phenotype"],
                  " .........."
                )
              )
              break
            }

            ## ----------------------------------------------------------
            ## Quandrant gate
            ## +/+
            ## ----------------------------------------------------------

            if (sign == "+,+") {
              ## define phenotpes
              fs[["phenotypes"]] <-
                ifelse(fs[[x]] > x_quant & fs[[y]] > y_quant,
                       master.table[mt, "phenotype"],
                       fs[["phenotypes"]])

              ## save gating results
              data.list[[toupper(master.table[mt, "createdGate"])]] <-
                fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

              if ((nrow(data.list[[toupper(master.table[mt, "createdGate"])]]) < 2)==FALSE) {

                ## percentage of gate
                per_gate <-
                  nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                  nrow(fs)

                ## labels
                # xLabel = x_quant + 3 * sd(fs[[x]])
                # yLabel = y_quant + 3 * sd(fs[[y]])
                xLabel = median(fs[fs[[x]] > x_quant & fs[[y]] > y_quant, ][[x]])
                yLabel = median(fs[fs[[x]] > x_quant & fs[[y]] > y_quant, ][[y]])

                ## plot
                p <- ggplot(fs, aes(x = .data[[x]],
                                    y = .data[[y]])) +
                  stat_density_2d(
                    aes(alpha = ..level..,
                        fill = ..level..),
                    size = 0,
                    bins = 50,
                    geom = "polygon"
                  ) +
                  xlab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(x))
                  )) +
                  ylab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(y))
                  )) +
                  geom_vline(
                    aes(xintercept = x_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_hline(
                    aes(yintercept = y_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = x_quant,
                      xmax = Inf,
                      ymin = y_quant,
                      ymax = Inf
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = xLabel,
                    y = yLabel,
                    label = paste0(round(per_gate, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  scale_fill_gradientn(colors = colfunc(200))  +
                  # annotation_logticks() +
                  geom_rug(size = 0.01) +
                  theme_bw() +
                  guides(alpha = "none",
                         fill = "none",
                         color = "none") +
                  theme(
                    axis.line = element_line(size = 0.75),
                    axis.text = element_text(size = 10,
                                             colour = "black"),
                    axis.title = element_text(size = 11,
                                              face = "italic")
                  ) +
                  theme(plot.title = element_text(face = "bold")) +
                  labs(title = paste0(master.table[mt, "phenotype"]),
                       subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"]))

                ## save plot
                facs.plot[[master.table[mt, "phenotype"]]] <- p

              }

            }

            ## ----------------------------------------------------------
            ## Quandrant gate
            ## -/+
            ## ----------------------------------------------------------

            if (sign == "-,+") {
              ## define phenotpes
              fs[["phenotypes"]] <-
                ifelse(fs[[x]] <= x_quant & fs[[y]] > y_quant,
                       master.table[mt, "phenotype"],
                       fs[["phenotypes"]])

              ## save gating results
              data.list[[toupper(master.table[mt, "createdGate"])]] <-
                fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

              if ((nrow(data.list[[toupper(master.table[mt, "createdGate"])]]) < 2)==FALSE) {

                ## percentage of gate
                per_gate <-
                  nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                  nrow(fs)

                ## label of gate
                # xLabel = x_quant + 3 * sd(fs[[x]])
                # yLabel = y_quant + 3 * sd(fs[[y]])
                xLabel = median(fs[fs[[x]] <= x_quant & fs[[y]] > y_quant, ][[x]])
                yLabel = median(fs[fs[[x]] <= x_quant & fs[[y]] > y_quant, ][[y]])

                # xLabel_per = x1+(x1/3)
                # yLabel_per = y1+(y1/3)
                ## plot
                p <- ggplot(fs, aes(x = .data[[x]],
                                    y = .data[[y]])) +
                  stat_density2d(
                    aes(alpha = ..level..,
                        fill = ..level..),
                    size = 0,
                    bins = 50,
                    geom = "polygon"
                  ) +
                  scale_fill_gradientn(colors = colfunc(200))  +
                  # annotation_logticks() +
                  guides(alpha = "none",
                         fill = "none",
                         color = "none") +
                  xlab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(x))
                  )) +
                  ylab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(y))
                  )) +
                  theme_bw() +
                  theme(
                    axis.line = element_line(size = 0.75),
                    axis.text = element_text(size = 10,
                                             colour = "black"),
                    axis.title = element_text(size = 11,
                                              face = "italic")
                  ) +
                  geom_vline(
                    aes(xintercept = x_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_hline(
                    aes(yintercept = y_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = -Inf,
                      xmax = x_quant,
                      ymin = y_quant,
                      ymax = Inf
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = xLabel,
                    y = yLabel,
                    label = paste0(round(per_gate, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  theme(plot.title = element_text(face = "bold")) +
                  labs(title = paste0(master.table[mt, "phenotype"]),
                       subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                  geom_rug(size = 0.01)

                ## save plot
                facs.plot[[master.table[mt, "phenotype"]]] <- p
              }

            }

            ## ----------------------------------------------------------
            ## Quandrant gate
            ## +/-
            ## ----------------------------------------------------------

            if (sign == "+,-") {
              ## define phenotpes
              fs[["phenotypes"]] <-
                ifelse(fs[[x]] > x_quant & fs[[y]] <= y_quant,
                       master.table[mt, "phenotype"],
                       fs[["phenotypes"]])

              ## save gating results
              data.list[[toupper(master.table[mt, "createdGate"])]] <-
                fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

              if ((nrow(data.list[[toupper(master.table[mt, "createdGate"])]]) < 2)==FALSE) {
                ## percentage of gate
                per_gate <-
                  nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                  nrow(fs)

                ## names of gate
                # xLabel = x_quant + 3 * sd(fs[[x]])
                # yLabel = y_quant + 3 * sd(fs[[y]])
                xLabel = median(fs[fs[[x]] > x_quant & fs[[y]] <= y_quant, ][[x]])
                yLabel = median(fs[fs[[x]] > x_quant & fs[[y]] <= y_quant, ][[y]])

                ## plot
                p <- ggplot(fs, aes(x = .data[[x]],
                                    y = .data[[y]])) +
                  stat_density2d(
                    aes(alpha = ..level..,
                        fill = ..level..),
                    size = 0,
                    bins = 50,
                    geom = "polygon"
                  ) +
                  scale_fill_gradientn(colors = colfunc(200))  +
                  # annotation_logticks() +
                  guides(alpha = "none",
                         fill = "none",
                         color = "none") +
                  xlab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(x))
                  )) +
                  ylab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(y))
                  )) +
                  theme_bw() +
                  theme(
                    axis.line = element_line(size = 0.75),
                    axis.text = element_text(size = 10,
                                             colour = "black"),
                    axis.title = element_text(size = 11,
                                              face = "italic")
                  ) +
                  geom_vline(
                    aes(xintercept = x_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_hline(
                    aes(yintercept = y_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = x_quant,
                      xmax = Inf,
                      ymin = -Inf,
                      ymax = y_quant
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = xLabel,
                    y = yLabel,
                    label = paste0(round(per_gate, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  theme(plot.title = element_text(face = "bold")) +
                  labs(title = paste0(master.table[mt, "phenotype"]),
                       subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                  geom_rug(size = 0.01)

                ## save plot
                facs.plot[[master.table[mt, "phenotype"]]] <- p
              }

            }
            ## ----------------------------------------------------------
            ## Quandrant gate
            ## -/-
            ## ----------------------------------------------------------
            if (sign == "-,-") {
              ## define phenotpes
              fs[["phenotypes"]] <-
                ifelse(fs[[x]] <= x_quant & fs[[y]] <= y_quant,
                       master.table[mt, "phenotype"],
                       fs[["phenotypes"]])

              ## save gating results
              data.list[[toupper(master.table[mt, "createdGate"])]] <-
                fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

              if ((nrow(data.list[[toupper(master.table[mt, "createdGate"])]]) < 3)==FALSE) {

                ## percentage of gate
                per_gate <-
                  nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                  nrow(fs)
                ## names of gate
                # xLabel = x_quant + 3 * sd(fs[[x]])
                # yLabel = y_quant + 3 * sd(fs[[y]])
                xLabel = median(fs[fs[[x]] <= x_quant & fs[[y]] <= y_quant, ][[x]])
                yLabel = median(fs[fs[[x]] <= x_quant & fs[[y]] <= y_quant, ][[y]])

                ## plot
                p <- ggplot(fs, aes(x = .data[[x]],
                                    y = .data[[y]])) +
                  stat_density2d(
                    aes(alpha = ..level..,
                        fill = ..level..),
                    size = 0,
                    bins = 50,
                    geom = "polygon"
                  ) +
                  # geom_point(alpha = 0.05,
                  #            size = 0.01) +
                  scale_fill_gradientn(colors = colfunc(200))  +
                  # annotation_logticks() +
                  guides(alpha = "none",
                         fill = "none",
                         color = "none") +
                  xlab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(x))
                  )) +
                  ylab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(y))
                  )) +
                  theme_bw() +
                  theme(
                    axis.line = element_line(size = 0.75),
                    axis.text = element_text(size = 10,
                                             colour = "black"),
                    axis.title = element_text(size = 11,
                                              face = "italic")
                  ) +
                  geom_vline(
                    aes(xintercept = x_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_hline(
                    aes(yintercept = y_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = -Inf,
                      xmax = x_quant,
                      ymin = -Inf,
                      ymax = y_quant
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = xLabel,
                    y = yLabel,
                    label = paste0(round(per_gate, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  theme(plot.title = element_text(face = "bold")) +
                  labs(title = paste0(master.table[mt, "phenotype"]),
                       subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                  geom_rug(size = 0.01)

                ## save plot
                facs.plot[[master.table[mt, "phenotype"]]] <- p
              }
            }

            ## ----------------------------------------------------------
            ## Quandrant gate
            ## */* for retaining all the gates
            ## ----------------------------------------------------------

            if (sign == "*,*") {
              for (q in 1:4) {
                ## ----------------------------------------------------------
                ## Quandrant gate
                ## -/- (Q1)
                ## ----------------------------------------------------------

                if (q == 1) {
                  ## set data frame
                  fs_q <- fs

                  ## define phenotpes
                  fs_q[["phenotypes"]] <-
                    ifelse(fs_q[[x]] <= x_quant &
                             fs_q[[y]] <= y_quant,
                           master.table[mt, "phenotype"],
                           fs_q[["phenotypes"]])

                  ## percentage of gate
                  per_gate_q1 <-
                    nrow(fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                    nrow(fs_q)

                  ## save gating results
                  data.list[[paste0(toupper(master.table[mt, "createdGate"]), "_", "Q", q)]] <-
                    fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]

                }

                ## ----------------------------------------------------------
                ## Quandrant gate
                ## -/+ (Q2)
                ## ----------------------------------------------------------
                else if (q == 2) {
                  ## set data frame
                  fs_q <- fs

                  ## define phenotpes
                  fs_q[["phenotypes"]] <-
                    ifelse(fs_q[[x]] <= x_quant &
                             fs_q[[y]] > y_quant,
                           master.table[mt, "phenotype"],
                           fs_q[["phenotypes"]])

                  ## percentage of gate
                  per_gate_q2 <-
                    nrow(fs[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                    nrow(fs_q)

                  ## save gating results
                  data.list[[paste0(toupper(master.table[mt, "createdGate"]), "_", "Q", q)]] <-
                    fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]
                }

                ## ----------------------------------------------------------
                ## Quandrant gate
                ## +/+ (Q3)
                ## ----------------------------------------------------------

                else if (q == 3) {
                  ## set data frame
                  fs_q <- fs

                  ## define phenotpes
                  fs_q[["phenotypes"]] <-
                    ifelse(fs_q[[x]] > x_quant &
                             fs_q[[y]] > y_quant,
                           master.table[mt, "phenotype"],
                           fs_q[["phenotypes"]])

                  ## percentage of gate
                  per_gate_q3 <-
                    nrow(fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                    nrow(fs_q)

                  ## save gating results
                  data.list[[paste0(toupper(master.table[mt, "createdGate"]), "_", "Q", q)]] <-
                    fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]
                }
                ## ----------------------------------------------------------
                ## Quandrant gate
                ## +/- (Q4)
                ## ----------------------------------------------------------
                else if (q == 4) {
                  ## set data frame
                  fs_q <- fs

                  ## define phenotpes
                  fs_q[["phenotypes"]] <-
                    ifelse(fs_q[[x]] > x_quant &
                             fs_q[[y]] <= y_quant,
                           master.table[mt, "phenotype"],
                           fs_q[["phenotypes"]])

                  ## percentage of gate
                  per_gate_q4 <-
                    nrow(fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                    nrow(fs_q)

                  ## save gating results
                  data.list[[paste0(toupper(master.table[mt, "createdGate"]), "_", "Q", q)]] <-
                    fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]
                }
              }
              if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

                ## plot
                p <- ggplot(fs, aes(x = .data[[x]],
                                    y = .data[[y]])) +
                  stat_density2d(
                    aes(alpha = ..level..,
                        fill = ..level..),
                    size = 0,
                    bins = 50,
                    geom = "polygon"
                  ) +
                  # geom_point(alpha = 0.05,
                  #            size = 0.01) +
                  scale_fill_gradientn(colors = colfunc(200))  +
                  # annotation_logticks() +
                  guides(alpha = "none",
                         fill = "none",
                         color = "none") +
                  xlab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(x))
                  )) +
                  ylab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(y))
                  )) +
                  theme_bw() +
                  theme(
                    axis.line = element_line(size = 0.75),
                    axis.text = element_text(size = 10,
                                             colour = "black"),
                    axis.title = element_text(size = 11,
                                              face = "italic")
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = -Inf,
                      xmax = x_quant,
                      ymin = -Inf,
                      ymax = y_quant
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = median(fs[fs[[x]] <= x_quant & fs[[y]] <= y_quant, ][[x]]),
                    y = median(fs[fs[[x]] <= x_quant & fs[[y]] <= y_quant, ][[y]]),
                    label = paste0(round(per_gate_q1, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = -Inf,
                      xmax = x_quant,
                      ymin = y_quant,
                      ymax = Inf
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = median(fs[fs[[x]] <= x_quant & fs[[y]] > y_quant, ][[x]]),
                    y = median(fs[fs[[x]] <= x_quant & fs[[y]] > y_quant, ][[y]]),
                    label = paste0(round(per_gate_q2, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = x_quant,
                      xmax = Inf,
                      ymin = y_quant,
                      ymax = Inf
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = median(fs[fs[[x]] > x_quant & fs[[y]] > y_quant, ][[x]]),
                    y = median(fs[fs[[x]] > x_quant & fs[[y]] > y_quant, ][[y]]),
                    label = paste0(round(per_gate_q3, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = x_quant,
                      xmax = Inf,
                      ymin = -Inf,
                      ymax = y_quant
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = median(fs[fs[[x]] > x_quant & fs[[y]] <= y_quant, ][[x]]),
                    y = median(fs[fs[[x]] > x_quant & fs[[y]] <= y_quant, ][[y]]),
                    label = paste0(round(per_gate_q4, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  theme(plot.title = element_text(face = "bold")) +
                  labs(title = paste0(master.table[mt, "phenotype"]),
                       subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                  geom_rug(size = 0.01)

                ## save plot
                facs.plot[[master.table[mt, "phenotype"]]] <- p
              }
            }
          }
          ## ----------------------------------------------------------
          ## rectangular gate
          ## No need to provide cutoff (default is 0.01 i.e 99.9% will be selected)
          ## If cutoff is provided, it will select 1-cutoff cells.
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "rectangleGate") {
            ## gate cutoff
            gateCutoff <- as.numeric(master.table[mt, "cutoff"])

            if (is.na(gateCutoff)) {
              gateCutoff <- 0.01
            }

            ## get quantile cutoff
            x_quant <-
              quantile(fs[[x]], c(gateCutoff, 1 - gateCutoff))
            y_quant <-
              quantile(fs[[y]], c(gateCutoff, 1 - gateCutoff))

            ## gate to choose
            x1 = x_quant[1]
            x2 = x_quant[2]
            y1 = y_quant[1]
            y2 = y_quant[2]

            ## define phenotpes
            fs[["phenotypes"]] <-
              ifelse(fs[[x]] > x1 &
                       fs[[y]] > y1 & fs[[x]] <= x2 & fs[[y]] <= y2,
                     master.table[mt, "phenotype"],
                     fs[["phenotypes"]])

            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              ## names of gate
              sign <-
                gsub(" ", "", master.table[mt, "markerDistribution"])
              char <- unlist(strsplit(sign, ","))
              name_gate <-
                paste0(x, "(", char[1], ")", y, "(", char[2], ")")
              
              xLabel = median(fs[[x]])
              yLabel = median(fs[[y]])

              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_rect(
                  mapping = aes(
                    xmin = x1,
                    xmax = x2,
                    ymin = y1,
                    ymax = y2
                  ),
                  color = "#3366ff",
                  fill = NA,
                  size = 0.25
                ) +
                annotate(
                  "label",
                  x = xLabel,
                  y = yLabel,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)

              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p

            }
          }

          ## ----------------------------------------------------------
          ## tail gate
          ## It is gate for one marker based on distribution....Only provide
          ## one markers, no need to provide cutoff
          ## ----------------------------------------------------------

          if (master.table[mt, "gateType"] %in% "tailGate") {
            if (length(xy) != 1) {
              print(
                paste0(
                  "please provide only 1 markers for ",
                  master.table[mt, "phenotype"],
                  " .........."
                )
              )
              break
            } else
              x <- xy

            ## gate cutoff
            gateCutoff <-
              mean(density(fs[[x]], kernel = "biweight")$x)

            ## get quantile cutoff
            x_quant <- quantile(fs[[x]], gateCutoff)

            ## names of gate
            sign <-
              gsub(" ", "", master.table[mt, "markerDistribution"])
            char <- unlist(strsplit(sign, ","))
            name_gate <- paste0(x, "(", char[1], ")")

            if (length(char) != 1) {
              print("please provide only 1 marker distribution")
              break
            }
            else if (char == "+") {
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] > x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              xLabel <- median(fs[fs[[x]] > x_quant, ][[x]])

              arrow.dir <- Inf
            }
            else if (char == "-") {
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] <= x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              xLabel <- median(fs[fs[[x]] <= x_quant, ][[x]])

              arrow.dir <- -Inf
            }

            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)

              ## plot
              p <- ggplot(fs, aes(x = .data[[x]])) +
                geom_density(fill = "grey", alpha = 0.5) +
                geom_segment(
                  aes(
                    x = x_quant,
                    y = 0,
                    xend = arrow.dir,
                    yend = 0
                  ),
                  lineend = "round",
                  linejoin = "round",
                  color = "#3366ff",
                  size = 0.25,
                  arrow = arrow(length = unit(0.2, "cm"))
                ) +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0("Proportion")) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_vline(aes(xintercept = x_quant),
                           size = 0.25,
                           color = "#3366ff") +
                annotate(
                  "label",
                  x = xLabel,
                  y = density(fs[[x]])$bw,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)

              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
            }

          }

          ## ----------------------------------------------------------
          ## quantile gate
          ## It is gate for one marker based on quantile distribution....Only provide
          ## one markers, do need to provide cutoff
          ## ----------------------------------------------------------

          if (master.table[mt, "gateType"] %in% "quantileGate") {
            if (length(xy) != 1) {
              print(
                paste0(
                  "please provide only 1 markers for ",
                  master.table[mt, "phenotype"],
                  " .........."
                )
              )
              break
            } else
              x <- xy

            ## gate cutoff
            gateCutoff <- as.numeric(master.table[mt, "cutoff"])

            if (is.na(gateCutoff)) {
              print("please provide quantile cutoff")
              break
            }

            ## get quantile cutoff
            x_quant <- gateCutoff

            ## names of gate
            sign <-
              gsub(" ", "", master.table[mt, "markerDistribution"])
            char <- unlist(strsplit(sign, ","))
            name_gate <- paste0(x, "(", char[1], ")")

            if (length(char) != 1) {
              print("please provide only 1 marker distribution")
              break
            } else if (char == "+") {
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] > x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              xLabel <- median(fs[fs[[x]] > x_quant, ][[x]])

              arrow.dir <- Inf

            } else if (char == "-") {
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] <= x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              xLabel <- median(fs[fs[[x]] <= x_quant, ][[x]])

              arrow.dir <- -Inf
            }

            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)

              ## plot
              p <- ggplot(fs, aes(x = .data[[x]])) +
                geom_density(fill = "grey", alpha = 0.5) +
                geom_segment(
                  aes(
                    x = x_quant,
                    y = 0,
                    xend = arrow.dir,
                    yend = 0
                  ),
                  lineend = "round",
                  linejoin = "round",
                  color = "#3366ff",
                  size = 0.25,
                  arrow = arrow(length = unit(0.1, "cm"))
                ) +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0("Proportion")) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_vline(aes(xintercept = x_quant),
                           size = 0.5,
                           color = "#3366ff") +
                annotate(
                  "label",
                  x = xLabel,
                  y = density(fs[[x]])$bw,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)

              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
            }
          }

          ## ----------------------------------------------------------
          ## hline gate
          ## It is gate for one marker based hline distribution....Only provide
          ## one markers and provide cutoff, ideal for Y variable
          ## ----------------------------------------------------------

          if (master.table[mt, "gateType"] %in% "hlineGate") {
            ## gate cutoff
            gateCutoff <- as.numeric(master.table[mt, "cutoff"])

            if (is.na(gateCutoff)) {
              print("pleasse provide cutoff")
              break
            }

            ## get quantile cutoff
            x_quant <- gateCutoff

            sign <-
              gsub(" ", "", master.table[mt, "markerDistribution"])
            char <- unlist(strsplit(sign, ","))
            name_gate <- paste0(x, "(", char[1], ")")

            if (length(char) != 1) {
              print("please provide the only one marker distribution")
            }
            if (char == "+") {
              ## gate to choose
              x1 = x_quant
              x2 = Inf
              y1 = -Inf
              y2 = Inf

              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] > x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)

              xLabel = median(fs[fs[[x]] > x_quant, ][[x]])
              yLabel = median(fs[fs[[x]] > x_quant, ][[y]])
            }
            if (char == "-") {
              ## gate to choose
              x1 = -Inf
              x2 = x_quant
              y1 = -Inf
              y2 = Inf

              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] <= x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)

              xLabel = median(fs[fs[[x]] <= x_quant, ][[x]])
              yLabel = median(fs[fs[[x]] <= x_quant, ][[y]])
            }

            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                # geom_rect(
                #   mapping = aes(
                #     xmin = x1,
                #     xmax = x2,
                #     ymin = y1,
                #     ymax = y2
                #   ),
                #   color = "#3366ff",
                #   fill = NA,
                #   size = 0.25
                # ) +
              geom_vline(xintercept = x_quant,color = "#3366ff",
                         fill = NA,
                         size = 0.5) +
                annotate(
                  "label",
                  x = xLabel,
                  y = yLabel,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)

              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p

            }
          }
          ## ----------------------------------------------------------
          ## vline gate
          ## It is gate for one marker based hline distribution....Only provide
          ## one markers and provide cutoff, ideal for X variable
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "vlineGate") {
            ## gate cutoff
            gateCutoff <- as.numeric(master.table[mt, "cutoff"])

            if (is.na(gateCutoff)) {
              print("pleasse provide cutoff")
              break
            }

            ## get quantile cutoff
            y_quant <- gateCutoff

            sign <-
              gsub(" ", "", master.table[mt, "markerDistribution"])
            char <- unlist(strsplit(sign, ","))
            name_gate <- paste0(y, "(", char[1], ")")

            if (length(char) != 1) {
              print("please provide the only one marker distribution")
            }
            else if (char == "+") {
              ## gate to choose
              x1 = -Inf
              x2 = Inf
              y1 = y_quant
              y2 = Inf

              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[y]] > y_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)

              xLabel = median(fs[fs[[y]] > y_quant, ][[x]])
              yLabel = median(fs[fs[[y]] > y_quant, ][[y]])
            }
            else if (char == "-") {
              ## gate to choose
              x1 = -Inf
              x2 = Inf
              y1 = -Inf
              y2 = y_quant

              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[y]] <= y_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              xLabel = median(fs[fs[[y]] <= y_quant, ][[x]])
              yLabel = median(fs[fs[[y]] <= y_quant, ][[y]])
            }

            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                # geom_rect(
                #   mapping = aes(
                #     xmin = x1,
                #     xmax = x2,
                #     ymin = y1,
                #     ymax = y2
                #   ),
                #   color = "#3366ff",
                #   fill = NA,
                #   size = 0.25
                # ) +
              geom_hline(yintercept = y_quant,color = "#3366ff",
                         fill = NA,
                         size = 0.5) +
                annotate(
                  "label",
                  x = xLabel,
                  y = yLabel,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)

              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p

            }

          }
          ## ----------------------------------------------------------
          ## boundary gate
          ## It is gate for selection of all the possible cells. No need to provide cutoffs
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "boundaryGate") {
            ## gate cutoff
            # minX <- median(fs[[x]]) - 3 * mad(fs[[x]], constant = 1)
            # maxX <- median(fs[[x]]) + 3 * mad(fs[[x]], constant = 1)
            #
            # minY <- median(fs[[y]]) - 3 * mad(fs[[y]], constant = 1)
            # maxY <- median(fs[[y]]) + 3 * mad(fs[[y]], constant = 1)
            minX <- min(fs[[x]])
            maxX <- max(fs[[x]])

            minY <- min(fs[[y]])
            maxY <- max(fs[[y]])

            ## gate to choose
            x1 = minX
            x2 = maxX
            y1 = minY
            y2 = maxY

            ## define phenotpes
            fs[["phenotypes"]] <-
              ifelse(fs[[x]] > x1 &
                       fs[[y]] > y1 & fs[[x]] <= x2 & fs[[y]] <= y2,
                     master.table[mt, "phenotype"],
                     fs[["phenotypes"]])

            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              ## names of gate
              sign <-
                gsub(" ", "", master.table[mt, "markerDistribution"])
              char <- unlist(strsplit(sign, ","))
              name_gate <-
                paste0(x, "(", char[1], ")", y, "(", char[2], ")")
              xLabel = median(fs[[x]])
              yLabel = median(fs[[y]])

              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_rect(
                  mapping = aes(
                    xmin = x1,
                    xmax = x2,
                    ymin = y1,
                    ymax = y2
                  ),
                  color = "#3366ff",
                  fill = NA,
                  size = 0.25
                ) +
                annotate(
                  "label",
                  x = xLabel,
                  y = yLabel,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)

              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p

            }

          }
          ## ----------------------------------------------------------
          ## singlet gate
          ## singlet gate fro selection of a singlet, no need to provide
          ## cutoff
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "singletGate") {
            formula <- as.formula(paste0(y, "~", x))
            model <- MASS::rlm(formula, fs, maxit = 5)

            which_min <- which.min(fs[[x]])
            which_max <- which.max(fs[[x]])
            x_extrema <-
              fs[c(which_min, which_max), colnames(fs) %in% c(x, y)]

            prediction_weights <- model$w[c(which_min, which_max)]

            predictions <-
              predict(model,
                      x_extrema,
                      interval = "prediction",
                      weights = prediction_weights)

            gate_vertices <-
              data.frame(rbind(
                cbind(x_extrema[[x]][1], predictions[1, "lwr"]),
                cbind(x_extrema[[x]][1], predictions[1, "upr"]),
                cbind(x_extrema[[x]][2], predictions[2, "upr"]),
                cbind(x_extrema[[x]][2], predictions[2, "lwr"])
              ))

            colnames(gate_vertices) <- c(x, y)

            ## get cells inside polygon
            polygon <- as.logical(sp::point.in.polygon(fs[[x]],
                                                   fs[[y]],
                                                   gate_vertices[[x]],
                                                   gate_vertices[[y]]))
            ## define phenotpes
            fs[["phenotypes"]] <- ifelse(polygon == TRUE,
                                         master.table[mt, "phenotype"],
                                         fs[["phenotypes"]])
            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              ## names of gate
              xLabel = mean(gate_vertices[[x]])
              yLabel = mean(gate_vertices[[y]])

              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_polygon(
                  data = gate_vertices,
                  aes(x = .data[[x]],
                      y = .data[[y]]),
                  color = "#3366ff",
                  fill = NA,
                  size = 0.25
                ) +
                annotate(
                  "label",
                  x = xLabel,
                  y = yLabel,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)

              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
            }
          }
          ## ----------------------------------------------------------
          ## # gate
          ## It will select whole poupation
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "#") {
            
            ## define phenotypes
              fs[["phenotypes"]] <- master.table[mt, "phenotype"]
              
            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
            
            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {
              
              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                )  +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)
              
              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
              
            }
            
          }
          ## ----------------------------------------------------------
          ## auto gate
          ## selection of a automated cluster, no need to provide
          ## cutoff...based on dbscan
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "autoGate") {
            datx <- fs[, colnames(fs) %in% c(x, y)]

            res = dbscan::dbscan(datx, eps = eps, minPts = nrow(fs) / 5)
            res = data.frame(datx, cluster = res$cluster)

            # Summarize each cluster
            info <- res %>%
              group_by(cluster) %>%
              summarise_all(., mean)

            info_coord <- info[info$cluster %in% K,]

            if (length(unique(res$cluster)) == 1) {
              print("no automated clusters found")
              break
            } else if (length(unique(res$cluster)) != 1) {
              for (c in unique(res$cluster)) {
                if (c == 0) {
                  res$cluster <- ifelse(res$cluster != 0, res$cluster, NA)
                } else if (c == K) {
                  ## define phenotpes
                  fs[["phenotypes"]] <- ifelse(res$cluster %in% c,
                                               master.table[mt, "phenotype"],
                                               fs[["phenotypes"]])
                  ## percentage of gate
                  per_gate <-
                    nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                    nrow(fs)

                  ## save gating results
                  data.list[[toupper(master.table[mt, "createdGate"])]] <-
                    fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
                }
              }
            }

            ## cluster
            data.sub <- res[res[["cluster"]] %in% K, ]
            kd <- ks::kde(data.sub[, c(x, y)], compute.cont = TRUE)
            res_cont <- with(
              kd,
              contourLines(
                x = eval.points[[1]],
                y = eval.points[[2]],
                z = estimate,
                levels = cont[scales::percent(5 /
                                                100)]
              )[[1]]
            )
            res_cont <- data.frame(res_cont)

            if ((is.null(res) || dim(ress)[1] < 2)==FALSE) {
            } else {
              ## plot
              p <- ggplot(res, aes(x = .data[[x]],
                                   y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_path(
                  data = res_cont,
                  aes(x = x,
                      y = y),
                  color = "#3366ff",
                  size = 0.25
                ) +
                annotate(
                  "label",
                  x = as.numeric(info_coord[1, x]),
                  y = as.numeric(info_coord[1, y]),
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)

              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
            }
          }
        }
      }
    }
    ## ----------------------------------------------------------
    ## other rows phenotypes
    ## ----------------------------------------------------------
    else {
      ## retrieve data for applied gate
      if (is.na(master.table[mt, "qoi"])) {
        fs <- data.list[[toupper(master.table[mt, "appliedGate"])]]
      } else {
        fs <-
          data.list[[paste0(toupper(master.table[mt, "appliedGate"]), "_", master.table[mt, "qoi"])]]
      }

      if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

        ### workaround for error in counter plots
        drop.cols <- c("IMAGENUMBER", "phenotypes")

        counter.errors <- fs %>%
          .[, setdiff(names(.), drop.cols)] %>%
          gather("variable", "value") %>%
          group_by(variable) %>%
          summarise(
            q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            q75 = quantile(value, probs = 0.75, na.rm = TRUE)
          ) %>%
          mutate(diff = q75 - q25) %>%
          filter(diff == 0)

        if (dim(counter.errors)[1] == 0) {
          ## select x,y
          xy <- toupper(master.table[mt, "markers"])
          xy <- gsub(" ", "", unlist(strsplit(xy, ",")))
          x <- xy[1]
          y <- xy[2]

          ## ----------------------------------------------------------
          ## Quandrant gate
          ## This gate is based on quandrant cutoff for X and Y (need
          ## to provide single cutoff for X and Y)
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "quadrantGate") {
            ## gate cutoff
            gateCutoff <- master.table[mt, "cutoff"]
            gateCutoff <-
              as.numeric(gsub(" ", "", unlist(strsplit(
                gateCutoff, ","
              ))))
            
            ##get quantile cutoff
            x_quant <- gateCutoff[1]
            y_quant <- gateCutoff[2]
            
            ## gate to choose
            sign <-
              gsub(" ", "", master.table[mt, "markerDistribution"])
            
            if (nchar(sign) != 3) {
              print(
                paste0(
                  "please provide binary marker distriubtion for ",
                  master.table[mt, "phenotype"],
                  " .........."
                )
              )
              break
            }
            
            ## ----------------------------------------------------------
            ## Quandrant gate
            ## +/+
            ## ----------------------------------------------------------
            
            if (sign == "+,+") {
              ## define phenotpes
              fs[["phenotypes"]] <-
                ifelse(fs[[x]] > x_quant & fs[[y]] > y_quant,
                       master.table[mt, "phenotype"],
                       fs[["phenotypes"]])
              
              ## save gating results
              data.list[[toupper(master.table[mt, "createdGate"])]] <-
                fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
              
              if ((nrow(data.list[[toupper(master.table[mt, "createdGate"])]]) < 2)==FALSE) {
                
                ## percentage of gate
                per_gate <-
                  nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                  nrow(fs)
                
                ## labels
                # xLabel = x_quant + 3 * sd(fs[[x]])
                # yLabel = y_quant + 3 * sd(fs[[y]])
                xLabel = median(fs[fs[[x]] > x_quant & fs[[y]] > y_quant, ][[x]])
                yLabel = median(fs[fs[[x]] > x_quant & fs[[y]] > y_quant, ][[y]])
                
                ## plot
                p <- ggplot(fs, aes(x = .data[[x]],
                                    y = .data[[y]])) +
                  stat_density_2d(
                    aes(alpha = ..level..,
                        fill = ..level..),
                    size = 0,
                    bins = 50,
                    geom = "polygon"
                  ) +
                  xlab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(x))
                  )) +
                  ylab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(y))
                  )) +
                  geom_vline(
                    aes(xintercept = x_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_hline(
                    aes(yintercept = y_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = x_quant,
                      xmax = Inf,
                      ymin = y_quant,
                      ymax = Inf
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = xLabel,
                    y = yLabel,
                    label = paste0(round(per_gate, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  scale_fill_gradientn(colors = colfunc(200))  +
                  # annotation_logticks() +
                  geom_rug(size = 0.01) +
                  theme_bw() +
                  guides(alpha = "none",
                         fill = "none",
                         color = "none") +
                  theme(
                    axis.line = element_line(size = 0.75),
                    axis.text = element_text(size = 10,
                                             colour = "black"),
                    axis.title = element_text(size = 11,
                                              face = "italic")
                  ) +
                  theme(plot.title = element_text(face = "bold")) +
                  labs(title = paste0(master.table[mt, "phenotype"]),
                       subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"]))
                
                ## save plot
                facs.plot[[master.table[mt, "phenotype"]]] <- p
                
              }
              
            }
            
            ## ----------------------------------------------------------
            ## Quandrant gate
            ## -/+
            ## ----------------------------------------------------------
            
            if (sign == "-,+") {
              ## define phenotpes
              fs[["phenotypes"]] <-
                ifelse(fs[[x]] <= x_quant & fs[[y]] > y_quant,
                       master.table[mt, "phenotype"],
                       fs[["phenotypes"]])
              
              ## save gating results
              data.list[[toupper(master.table[mt, "createdGate"])]] <-
                fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
              
              if ((nrow(data.list[[toupper(master.table[mt, "createdGate"])]]) < 2)==FALSE) {
                
                ## percentage of gate
                per_gate <-
                  nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                  nrow(fs)
                
                ## label of gate
                # xLabel = x_quant + 3 * sd(fs[[x]])
                # yLabel = y_quant + 3 * sd(fs[[y]])
                xLabel = median(fs[fs[[x]] <= x_quant & fs[[y]] > y_quant, ][[x]])
                yLabel = median(fs[fs[[x]] <= x_quant & fs[[y]] > y_quant, ][[y]])
                
                # xLabel_per = x1+(x1/3)
                # yLabel_per = y1+(y1/3)
                ## plot
                p <- ggplot(fs, aes(x = .data[[x]],
                                    y = .data[[y]])) +
                  stat_density2d(
                    aes(alpha = ..level..,
                        fill = ..level..),
                    size = 0,
                    bins = 50,
                    geom = "polygon"
                  ) +
                  scale_fill_gradientn(colors = colfunc(200))  +
                  # annotation_logticks() +
                  guides(alpha = "none",
                         fill = "none",
                         color = "none") +
                  xlab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(x))
                  )) +
                  ylab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(y))
                  )) +
                  theme_bw() +
                  theme(
                    axis.line = element_line(size = 0.75),
                    axis.text = element_text(size = 10,
                                             colour = "black"),
                    axis.title = element_text(size = 11,
                                              face = "italic")
                  ) +
                  geom_vline(
                    aes(xintercept = x_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_hline(
                    aes(yintercept = y_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = -Inf,
                      xmax = x_quant,
                      ymin = y_quant,
                      ymax = Inf
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = xLabel,
                    y = yLabel,
                    label = paste0(round(per_gate, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  theme(plot.title = element_text(face = "bold")) +
                  labs(title = paste0(master.table[mt, "phenotype"]),
                       subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                  geom_rug(size = 0.01)
                
                ## save plot
                facs.plot[[master.table[mt, "phenotype"]]] <- p
              }
              
            }
            
            ## ----------------------------------------------------------
            ## Quandrant gate
            ## +/-
            ## ----------------------------------------------------------
            
            if (sign == "+,-") {
              ## define phenotpes
              fs[["phenotypes"]] <-
                ifelse(fs[[x]] > x_quant & fs[[y]] <= y_quant,
                       master.table[mt, "phenotype"],
                       fs[["phenotypes"]])
              
              ## save gating results
              data.list[[toupper(master.table[mt, "createdGate"])]] <-
                fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
              
              if ((nrow(data.list[[toupper(master.table[mt, "createdGate"])]]) < 2)==FALSE) {
                ## percentage of gate
                per_gate <-
                  nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                  nrow(fs)
                
                ## names of gate
                # xLabel = x_quant + 3 * sd(fs[[x]])
                # yLabel = y_quant + 3 * sd(fs[[y]])
                xLabel = median(fs[fs[[x]] > x_quant & fs[[y]] <= y_quant, ][[x]])
                yLabel = median(fs[fs[[x]] > x_quant & fs[[y]] <= y_quant, ][[y]])
                
                ## plot
                p <- ggplot(fs, aes(x = .data[[x]],
                                    y = .data[[y]])) +
                  stat_density2d(
                    aes(alpha = ..level..,
                        fill = ..level..),
                    size = 0,
                    bins = 50,
                    geom = "polygon"
                  ) +
                  scale_fill_gradientn(colors = colfunc(200))  +
                  # annotation_logticks() +
                  guides(alpha = "none",
                         fill = "none",
                         color = "none") +
                  xlab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(x))
                  )) +
                  ylab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(y))
                  )) +
                  theme_bw() +
                  theme(
                    axis.line = element_line(size = 0.75),
                    axis.text = element_text(size = 10,
                                             colour = "black"),
                    axis.title = element_text(size = 11,
                                              face = "italic")
                  ) +
                  geom_vline(
                    aes(xintercept = x_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_hline(
                    aes(yintercept = y_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = x_quant,
                      xmax = Inf,
                      ymin = -Inf,
                      ymax = y_quant
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = xLabel,
                    y = yLabel,
                    label = paste0(round(per_gate, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  theme(plot.title = element_text(face = "bold")) +
                  labs(title = paste0(master.table[mt, "phenotype"]),
                       subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                  geom_rug(size = 0.01)
                
                ## save plot
                facs.plot[[master.table[mt, "phenotype"]]] <- p
              }
              
            }
            ## ----------------------------------------------------------
            ## Quandrant gate
            ## -/-
            ## ----------------------------------------------------------
            if (sign == "-,-") {
              ## define phenotpes
              fs[["phenotypes"]] <-
                ifelse(fs[[x]] <= x_quant & fs[[y]] <= y_quant,
                       master.table[mt, "phenotype"],
                       fs[["phenotypes"]])
              
              ## save gating results
              data.list[[toupper(master.table[mt, "createdGate"])]] <-
                fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
              
              if ((nrow(data.list[[toupper(master.table[mt, "createdGate"])]]) < 3)==FALSE) {
                
                ## percentage of gate
                per_gate <-
                  nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                  nrow(fs)
                ## names of gate
                # xLabel = x_quant + 3 * sd(fs[[x]])
                # yLabel = y_quant + 3 * sd(fs[[y]])
                xLabel = median(fs[fs[[x]] <= x_quant & fs[[y]] <= y_quant, ][[x]])
                yLabel = median(fs[fs[[x]] <= x_quant & fs[[y]] <= y_quant, ][[y]])
                
                ## plot
                p <- ggplot(fs, aes(x = .data[[x]],
                                    y = .data[[y]])) +
                  stat_density2d(
                    aes(alpha = ..level..,
                        fill = ..level..),
                    size = 0,
                    bins = 50,
                    geom = "polygon"
                  ) +
                  # geom_point(alpha = 0.05,
                  #            size = 0.01) +
                  scale_fill_gradientn(colors = colfunc(200))  +
                  # annotation_logticks() +
                  guides(alpha = "none",
                         fill = "none",
                         color = "none") +
                  xlab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(x))
                  )) +
                  ylab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(y))
                  )) +
                  theme_bw() +
                  theme(
                    axis.line = element_line(size = 0.75),
                    axis.text = element_text(size = 10,
                                             colour = "black"),
                    axis.title = element_text(size = 11,
                                              face = "italic")
                  ) +
                  geom_vline(
                    aes(xintercept = x_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_hline(
                    aes(yintercept = y_quant),
                    color = "#377eb8",
                    size = 0.5,
                    lty = "dotted"
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = -Inf,
                      xmax = x_quant,
                      ymin = -Inf,
                      ymax = y_quant
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = xLabel,
                    y = yLabel,
                    label = paste0(round(per_gate, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  theme(plot.title = element_text(face = "bold")) +
                  labs(title = paste0(master.table[mt, "phenotype"]),
                       subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                  geom_rug(size = 0.01)
                
                ## save plot
                facs.plot[[master.table[mt, "phenotype"]]] <- p
              }
            }
            
            ## ----------------------------------------------------------
            ## Quandrant gate
            ## */* for retaining all the gates
            ## ----------------------------------------------------------
            
            if (sign == "*,*") {
              for (q in 1:4) {
                ## ----------------------------------------------------------
                ## Quandrant gate
                ## -/- (Q1)
                ## ----------------------------------------------------------
                
                if (q == 1) {
                  ## set data frame
                  fs_q <- fs
                  
                  ## define phenotpes
                  fs_q[["phenotypes"]] <-
                    ifelse(fs_q[[x]] <= x_quant &
                             fs_q[[y]] <= y_quant,
                           master.table[mt, "phenotype"],
                           fs_q[["phenotypes"]])
                  
                  ## percentage of gate
                  per_gate_q1 <-
                    nrow(fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                    nrow(fs_q)
                  
                  ## save gating results
                  data.list[[paste0(toupper(master.table[mt, "createdGate"]), "_", "Q", q)]] <-
                    fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]
                  
                }
                
                ## ----------------------------------------------------------
                ## Quandrant gate
                ## -/+ (Q2)
                ## ----------------------------------------------------------
                else if (q == 2) {
                  ## set data frame
                  fs_q <- fs
                  
                  ## define phenotpes
                  fs_q[["phenotypes"]] <-
                    ifelse(fs_q[[x]] <= x_quant &
                             fs_q[[y]] > y_quant,
                           master.table[mt, "phenotype"],
                           fs_q[["phenotypes"]])
                  
                  ## percentage of gate
                  per_gate_q2 <-
                    nrow(fs[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                    nrow(fs_q)
                  
                  ## save gating results
                  data.list[[paste0(toupper(master.table[mt, "createdGate"]), "_", "Q", q)]] <-
                    fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]
                }
                
                ## ----------------------------------------------------------
                ## Quandrant gate
                ## +/+ (Q3)
                ## ----------------------------------------------------------
                
                else if (q == 3) {
                  ## set data frame
                  fs_q <- fs
                  
                  ## define phenotpes
                  fs_q[["phenotypes"]] <-
                    ifelse(fs_q[[x]] > x_quant &
                             fs_q[[y]] > y_quant,
                           master.table[mt, "phenotype"],
                           fs_q[["phenotypes"]])
                  
                  ## percentage of gate
                  per_gate_q3 <-
                    nrow(fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                    nrow(fs_q)
                  
                  ## save gating results
                  data.list[[paste0(toupper(master.table[mt, "createdGate"]), "_", "Q", q)]] <-
                    fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]
                }
                ## ----------------------------------------------------------
                ## Quandrant gate
                ## +/- (Q4)
                ## ----------------------------------------------------------
                else if (q == 4) {
                  ## set data frame
                  fs_q <- fs
                  
                  ## define phenotpes
                  fs_q[["phenotypes"]] <-
                    ifelse(fs_q[[x]] > x_quant &
                             fs_q[[y]] <= y_quant,
                           master.table[mt, "phenotype"],
                           fs_q[["phenotypes"]])
                  
                  ## percentage of gate
                  per_gate_q4 <-
                    nrow(fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                    nrow(fs_q)
                  
                  ## save gating results
                  data.list[[paste0(toupper(master.table[mt, "createdGate"]), "_", "Q", q)]] <-
                    fs_q[fs_q[["phenotypes"]] %in% master.table[mt, "phenotype"],]
                }
              }
              if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {
                
                ## plot
                p <- ggplot(fs, aes(x = .data[[x]],
                                    y = .data[[y]])) +
                  stat_density2d(
                    aes(alpha = ..level..,
                        fill = ..level..),
                    size = 0,
                    bins = 50,
                    geom = "polygon"
                  ) +
                  # geom_point(alpha = 0.05,
                  #            size = 0.01) +
                  scale_fill_gradientn(colors = colfunc(200))  +
                  # annotation_logticks() +
                  guides(alpha = "none",
                         fill = "none",
                         color = "none") +
                  xlab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(x))
                  )) +
                  ylab(paste0(
                    "scaled mean intensity \n",
                    tools::toTitleCase(tolower(y))
                  )) +
                  theme_bw() +
                  theme(
                    axis.line = element_line(size = 0.75),
                    axis.text = element_text(size = 10,
                                             colour = "black"),
                    axis.title = element_text(size = 11,
                                              face = "italic")
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = -Inf,
                      xmax = x_quant,
                      ymin = -Inf,
                      ymax = y_quant
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = median(fs[fs[[x]] <= x_quant & fs[[y]] <= y_quant, ][[x]]),
                    y = median(fs[fs[[x]] <= x_quant & fs[[y]] <= y_quant, ][[y]]),
                    label = paste0(round(per_gate_q1, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = -Inf,
                      xmax = x_quant,
                      ymin = y_quant,
                      ymax = Inf
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = median(fs[fs[[x]] <= x_quant & fs[[y]] > y_quant, ][[x]]),
                    y = median(fs[fs[[x]] <= x_quant & fs[[y]] > y_quant, ][[y]]),
                    label = paste0(round(per_gate_q2, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = x_quant,
                      xmax = Inf,
                      ymin = y_quant,
                      ymax = Inf
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = median(fs[fs[[x]] > x_quant & fs[[y]] > y_quant, ][[x]]),
                    y = median(fs[fs[[x]] > x_quant & fs[[y]] > y_quant, ][[y]]),
                    label = paste0(round(per_gate_q3, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  geom_rect(
                    mapping = aes(
                      xmin = x_quant,
                      xmax = Inf,
                      ymin = -Inf,
                      ymax = y_quant
                    ),
                    color = "#3366ff",
                    fill = NA,
                    size = 0.5
                  ) +
                  annotate(
                    "label",
                    x = median(fs[fs[[x]] > x_quant & fs[[y]] <= y_quant, ][[x]]),
                    y = median(fs[fs[[x]] > x_quant & fs[[y]] <= y_quant, ][[y]]),
                    label = paste0(round(per_gate_q4, 2), "%"),
                    fontface = 2,
                    size = 3,
                    fill=NA
                  ) +
                  theme(plot.title = element_text(face = "bold")) +
                  labs(title = paste0(master.table[mt, "phenotype"]),
                       subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                  geom_rug(size = 0.01)
                
                ## save plot
                facs.plot[[master.table[mt, "phenotype"]]] <- p
              }
            }
          }
          ## ----------------------------------------------------------
          ## rectangular gate
          ## No need to provide cutoff (default is 0.01 i.e 99.9% will be selected)
          ## If cutoff is provided, it will select 1-cutoff cells.
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "rectangleGate") {
            ## gate cutoff
            gateCutoff <- as.numeric(master.table[mt, "cutoff"])
            
            if (is.na(gateCutoff)) {
              gateCutoff <- 0.01
            }
            
            ## get quantile cutoff
            x_quant <-
              quantile(fs[[x]], c(gateCutoff, 1 - gateCutoff))
            y_quant <-
              quantile(fs[[y]], c(gateCutoff, 1 - gateCutoff))
            
            ## gate to choose
            x1 = x_quant[1]
            x2 = x_quant[2]
            y1 = y_quant[1]
            y2 = y_quant[2]
            
            ## define phenotpes
            fs[["phenotypes"]] <-
              ifelse(fs[[x]] > x1 &
                       fs[[y]] > y1 & fs[[x]] <= x2 & fs[[y]] <= y2,
                     master.table[mt, "phenotype"],
                     fs[["phenotypes"]])
            
            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
            
            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {
              
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              ## names of gate
              sign <-
                gsub(" ", "", master.table[mt, "markerDistribution"])
              char <- unlist(strsplit(sign, ","))
              name_gate <-
                paste0(x, "(", char[1], ")", y, "(", char[2], ")")
              
              xLabel = median(fs[[x]])
              yLabel = median(fs[[y]])
              
              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_rect(
                  mapping = aes(
                    xmin = x1,
                    xmax = x2,
                    ymin = y1,
                    ymax = y2
                  ),
                  color = "#3366ff",
                  fill = NA,
                  size = 0.25
                ) +
                annotate(
                  "label",
                  x = xLabel,
                  y = yLabel,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)
              
              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
              
            }
          }
          ## ----------------------------------------------------------
          ## tail gate
          ## It is gate for one marker based on distribution....Only provide
          ## one markers, no need to provide cutoff
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "tailGate") {
            if (length(xy) != 1) {
              print(
                paste0(
                  "please provide only 1 markers for ",
                  master.table[mt, "phenotype"],
                  " .........."
                )
              )
              break
            } else
              x <- xy
            
            ## gate cutoff
            gateCutoff <-
              mean(density(fs[[x]], kernel = "biweight")$x)
            
            ## get quantile cutoff
            x_quant <- quantile(fs[[x]], gateCutoff)
            
            ## names of gate
            sign <-
              gsub(" ", "", master.table[mt, "markerDistribution"])
            char <- unlist(strsplit(sign, ","))
            name_gate <- paste0(x, "(", char[1], ")")
            
            if (length(char) != 1) {
              print("please provide only 1 marker distribution")
              break
            }
            else if (char == "+") {
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] > x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              xLabel <- median(fs[fs[[x]] > x_quant, ][[x]])
              
              arrow.dir <- Inf
            }
            else if (char == "-") {
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] <= x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              xLabel <- median(fs[fs[[x]] <= x_quant, ][[x]])
              
              arrow.dir <- -Inf
            }
            
            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
            
            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {
              
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              
              ## plot
              p <- ggplot(fs, aes(x = .data[[x]])) +
                geom_density(fill = "grey", alpha = 0.5) +
                geom_segment(
                  aes(
                    x = x_quant,
                    y = 0,
                    xend = arrow.dir,
                    yend = 0
                  ),
                  lineend = "round",
                  linejoin = "round",
                  color = "#3366ff",
                  size = 0.25,
                  arrow = arrow(length = unit(0.2, "cm"))
                ) +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0("Proportion")) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_vline(aes(xintercept = x_quant),
                           size = 0.25,
                           color = "#3366ff") +
                annotate(
                  "label",
                  x = xLabel,
                  y = density(fs[[x]])$bw,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)
              
              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
            }
            
          }
          ## ----------------------------------------------------------
          ## quantile gate
          ## It is gate for one marker based on quantile distribution....Only provide
          ## one markers, do need to provide cutoff
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "quantileGate") {
            if (length(xy) != 1) {
              print(
                paste0(
                  "please provide only 1 markers for ",
                  master.table[mt, "phenotype"],
                  " .........."
                )
              )
              break
            } else
              x <- xy
            
            ## gate cutoff
            gateCutoff <- as.numeric(master.table[mt, "cutoff"])
            
            if (is.na(gateCutoff)) {
              print("please provide quantile cutoff")
              break
            }
            
            ## get quantile cutoff
            x_quant <- gateCutoff
            
            ## names of gate
            sign <-
              gsub(" ", "", master.table[mt, "markerDistribution"])
            char <- unlist(strsplit(sign, ","))
            name_gate <- paste0(x, "(", char[1], ")")
            
            if (length(char) != 1) {
              print("please provide only 1 marker distribution")
              break
            } else if (char == "+") {
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] > x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              xLabel <- median(fs[fs[[x]] > x_quant, ][[x]])
              
              arrow.dir <- Inf
              
            } else if (char == "-") {
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] <= x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              xLabel <- median(fs[fs[[x]] <= x_quant, ][[x]])
              
              arrow.dir <- -Inf
            }
            
            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
            
            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {
              
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              
              ## plot
              p <- ggplot(fs, aes(x = .data[[x]])) +
                geom_density(fill = "grey", alpha = 0.5) +
                geom_segment(
                  aes(
                    x = x_quant,
                    y = 0,
                    xend = arrow.dir,
                    yend = 0
                  ),
                  lineend = "round",
                  linejoin = "round",
                  color = "#3366ff",
                  size = 0.25,
                  arrow = arrow(length = unit(0.1, "cm"))
                ) +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0("Proportion")) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_vline(aes(xintercept = x_quant),
                           size = 0.5,
                           color = "#3366ff") +
                annotate(
                  "label",
                  x = xLabel,
                  y = density(fs[[x]])$bw,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)
              
              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
            }
          }
          ## ----------------------------------------------------------
          ## hline gate
          ## It is gate for one marker based hline distribution....Only provide
          ## one markers and provide cutoff, ideal for X variable
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "hlineGate") {
            ## gate cutoff
            gateCutoff <- as.numeric(master.table[mt, "cutoff"])
            
            if (is.na(gateCutoff)) {
              print("pleasse provide cutoff")
              break
            }
            
            ## get quantile cutoff
            x_quant <- gateCutoff
            
            sign <-
              gsub(" ", "", master.table[mt, "markerDistribution"])
            char <- unlist(strsplit(sign, ","))
            name_gate <- paste0(x, "(", char[1], ")")
            
            if (length(char) != 1) {
              print("please provide the only one marker distribution")
            }
            if (char == "+") {
              ## gate to choose
              x1 = x_quant
              x2 = Inf
              y1 = -Inf
              y2 = Inf
              
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] > x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              
              xLabel = median(fs[fs[[x]] > x_quant, ][[x]])
              yLabel = median(fs[fs[[x]] > x_quant, ][[y]])
            }
            if (char == "-") {
              ## gate to choose
              x1 = -Inf
              x2 = x_quant
              y1 = -Inf
              y2 = Inf
              
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[x]] <= x_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              
              xLabel = median(fs[fs[[x]] <= x_quant, ][[x]])
              yLabel = median(fs[fs[[x]] <= x_quant, ][[y]])
            }
            
            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
            
            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {
              
              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                # geom_rect(
                #   mapping = aes(
                #     xmin = x1,
                #     xmax = x2,
                #     ymin = y1,
                #     ymax = y2
                #   ),
                #   color = "#3366ff",
                #   fill = NA,
                #   size = 0.25
                # ) +
              geom_vline(xintercept = x_quant,color = "#3366ff",
                         fill = NA,
                         size = 0.5) +
                annotate(
                  "label",
                  x = xLabel,
                  y = yLabel,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)
              
              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
              
            }
          }
          ## ----------------------------------------------------------
          ## vline gate
          ## It is gate for one marker based hline distribution....Only provide
          ## one markers and provide cutoff, ideal for Y variable
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "vlineGate") {
            ## gate cutoff
            gateCutoff <- as.numeric(master.table[mt, "cutoff"])
            
            if (is.na(gateCutoff)) {
              print("pleasse provide cutoff")
              break
            }
            
            ## get quantile cutoff
            y_quant <- gateCutoff
            
            sign <-
              gsub(" ", "", master.table[mt, "markerDistribution"])
            char <- unlist(strsplit(sign, ","))
            name_gate <- paste0(y, "(", char[1], ")")
            
            if (length(char) != 1) {
              print("please provide the only one marker distribution")
            }
            else if (char == "+") {
              ## gate to choose
              x1 = -Inf
              x2 = Inf
              y1 = y_quant
              y2 = Inf
              
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[y]] > y_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              
              xLabel = median(fs[fs[[y]] > y_quant, ][[x]])
              yLabel = median(fs[fs[[y]] > y_quant, ][[y]])
            }
            else if (char == "-") {
              ## gate to choose
              x1 = -Inf
              x2 = Inf
              y1 = -Inf
              y2 = y_quant
              
              ## define phenotpes
              fs[["phenotypes"]] <- ifelse(fs[[y]] <= y_quant,
                                           master.table[mt, "phenotype"],
                                           fs[["phenotypes"]])
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              xLabel = median(fs[fs[[y]] <= y_quant, ][[x]])
              yLabel = median(fs[fs[[y]] <= y_quant, ][[y]])
            }
            
            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
            
            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {
              
              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                # geom_rect(
                #   mapping = aes(
                #     xmin = x1,
                #     xmax = x2,
                #     ymin = y1,
                #     ymax = y2
                #   ),
                #   color = "#3366ff",
                #   fill = NA,
                #   size = 0.25
                # ) +
              geom_hline(yintercept = y_quant,color = "#3366ff",
                         fill = NA,
                         size = 0.5) +
                annotate(
                  "label",
                  x = xLabel,
                  y = yLabel,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)
              
              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
              
            }
            
          }
          ## ----------------------------------------------------------
          ## boundary gate
          ## It is gate for selection of all the possible cells. No need to provide cutoffs
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "boundaryGate") {
            ## gate cutoff
            # minX <- median(fs[[x]]) - 3 * mad(fs[[x]], constant = 1)
            # maxX <- median(fs[[x]]) + 3 * mad(fs[[x]], constant = 1)
            #
            # minY <- median(fs[[y]]) - 3 * mad(fs[[y]], constant = 1)
            # maxY <- median(fs[[y]]) + 3 * mad(fs[[y]], constant = 1)
            minX <- min(fs[[x]])
            maxX <- max(fs[[x]])
            
            minY <- min(fs[[y]])
            maxY <- max(fs[[y]])
            
            ## gate to choose
            x1 = minX
            x2 = maxX
            y1 = minY
            y2 = maxY
            
            ## define phenotpes
            fs[["phenotypes"]] <-
              ifelse(fs[[x]] > x1 &
                       fs[[y]] > y1 & fs[[x]] <= x2 & fs[[y]] <= y2,
                     master.table[mt, "phenotype"],
                     fs[["phenotypes"]])
            
            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
            
            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {
              
              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              ## names of gate
              sign <-
                gsub(" ", "", master.table[mt, "markerDistribution"])
              char <- unlist(strsplit(sign, ","))
              name_gate <-
                paste0(x, "(", char[1], ")", y, "(", char[2], ")")
              xLabel = median(fs[[x]])
              yLabel = median(fs[[y]])
              
              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_rect(
                  mapping = aes(
                    xmin = x1,
                    xmax = x2,
                    ymin = y1,
                    ymax = y2
                  ),
                  color = "#3366ff",
                  fill = NA,
                  size = 0.25
                ) +
                annotate(
                  "label",
                  x = xLabel,
                  y = yLabel,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)
              
              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
              
            }
            
          }
          ## ----------------------------------------------------------
          ## singlet gate
          ## singlet gate fro selection of a singlet, no need to provide
          ## cutoff
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "singletGate") {
            formula <- as.formula(paste0(y, "~", x))
            model <- MASS::rlm(formula, fs, maxit = 5)

            which_min <- which.min(fs[[x]])
            which_max <- which.max(fs[[x]])
            x_extrema <-
              fs[c(which_min, which_max), colnames(fs) %in% c(x, y)]

            prediction_weights <- model$w[c(which_min, which_max)]

            predictions <-
              predict(model,
                      x_extrema,
                      interval = "prediction",
                      weights = prediction_weights)

            gate_vertices <-
              data.frame(rbind(
                cbind(x_extrema[[x]][1], predictions[1, "lwr"]),
                cbind(x_extrema[[x]][1], predictions[1, "upr"]),
                cbind(x_extrema[[x]][2], predictions[2, "upr"]),
                cbind(x_extrema[[x]][2], predictions[2, "lwr"])
              ))

            colnames(gate_vertices) <- c(x, y)

            ## get cells inside polygon
            polygon <- as.logical(sp::point.in.polygon(fs[[x]],
                                                   fs[[y]],
                                                   gate_vertices[[x]],
                                                   gate_vertices[[y]]))
            ## define phenotpes
            fs[["phenotypes"]] <- ifelse(polygon == TRUE,
                                         master.table[mt, "phenotype"],
                                         fs[["phenotypes"]])

            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]

            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {

              ## percentage of gate
              per_gate <-
                nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                nrow(fs)
              ## names of gate
              xLabel = mean(gate_vertices[[x]])
              yLabel = mean(gate_vertices[[y]])

              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_polygon(
                  data = gate_vertices,
                  aes(x = .data[[x]],
                      y = .data[[y]]),
                  color = "#3366ff",
                  fill = NA,
                  size = 0.25
                ) +
                annotate(
                  "label",
                  x = xLabel,
                  y = yLabel,
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)

              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p

            }

          }
          ## ----------------------------------------------------------
          ## # gate
          ## It will select whole population
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "#") {
            
            ## define phenotypes
            fs[["phenotypes"]] <- master.table[mt, "phenotype"]
            
            ## save gating results
            data.list[[toupper(master.table[mt, "createdGate"])]] <-
              fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
            
            if ((is.null(fs) || dim(fs)[1] < 2)==FALSE) {
              
              ## plot
              p <- ggplot(fs, aes(x = .data[[x]],
                                  y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                )  +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)
              
              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
              
            }
            
          }
          ## ----------------------------------------------------------
          ## auto gate
          ## selection of a automated cluster, no need to provide
          ## cutoff...based on dbscan
          ## ----------------------------------------------------------
          if (master.table[mt, "gateType"] %in% "autoGate") {
            datx <- fs[, colnames(fs) %in% c(x, y)]

            res = dbscan::dbscan(datx, eps = eps, minPts = nrow(fs) / 10)
            res = data.frame(datx, cluster = res$cluster)

            # Summarize each cluster
            info <- res %>%
              group_by(cluster) %>%
              summarise_all(., mean)

            info_coord <- info[info$cluster %in% K,]

            if (length(unique(res$cluster)) == 1) {
              print("no automated clusters found")
              break
            } else if (length(unique(res$cluster)) != 1) {
              for (c in unique(res$cluster)) {
                if (c == 0) {
                  res$cluster <- ifelse(res$cluster != 0, res$cluster, NA)
                } else if (c == K) {
                  ## define phenotpes
                  fs[["phenotypes"]] <- ifelse(res$cluster %in% c,
                                               master.table[mt, "phenotype"],
                                               fs[["phenotypes"]])
                  ## percentage of gate
                  per_gate <-
                    nrow(fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]) * 100 /
                    nrow(fs)

                  ## save gating results
                  data.list[[toupper(master.table[mt, "createdGate"])]] <-
                    fs[fs[["phenotypes"]] %in% master.table[mt, "phenotype"],]
                }
              }
            }

            ## cluster
            data.sub <- res[res[["cluster"]] %in% K, ]
            kd <- ks::kde(data.sub[, c(x, y)], compute.cont = TRUE)
            res_cont <- with(
              kd,
              contourLines(
                x = eval.points[[1]],
                y = eval.points[[2]],
                z = estimate,
                levels = cont[scales::percent(5 /
                                                100)]
              )[[1]]
            )
            res_cont <- data.frame(res_cont)

            if ((is.null(res) || dim(res)[1] < 2)==FALSE) {

              ## plot
              p <- ggplot(res, aes(x = .data[[x]],
                                   y = .data[[y]])) +
                stat_density2d(
                  aes(alpha = ..level..,
                      fill = ..level..),
                  size = 0,
                  bins = 50,
                  geom = "polygon"
                ) +
                # geom_point(alpha = 0.05,
                #            size = 0.01) +
                scale_fill_gradientn(colors = colfunc(200))  +
                # annotation_logticks() +
                guides(alpha = "none",
                       fill = "none",
                       color = "none") +
                xlab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(x))
                )) +
                ylab(paste0(
                  "scaled mean intensity \n",
                  tools::toTitleCase(tolower(y))
                )) +
                theme_bw() +
                theme(
                  axis.line = element_line(size = 0.75),
                  axis.text = element_text(size = 10,
                                           colour = "black"),
                  axis.title = element_text(size = 11,
                                            face = "italic")
                ) +
                geom_path(
                  data = res_cont,
                  aes(x = x,
                      y = y),
                  color = "#3366ff",
                  size = 0.25
                ) +
                annotate(
                  "label",
                  x = as.numeric(info_coord[1, x]),
                  y = as.numeric(info_coord[1, y]),
                  label = paste0(round(per_gate, 2), "%"),
                  fontface = 2,
                  size = 3,
                  fill=NA
                ) +
                theme(plot.title = element_text(face = "bold")) +
                labs(title = paste0(master.table[mt, "phenotype"]),
                     subtitle = paste(master.table[mt, "appliedGate"], master.table[mt, "qoi"])) +
                geom_rug(size = 0.01)

              ## save plot
              facs.plot[[master.table[mt, "phenotype"]]] <- p
            }
          }
        }
      }
    }
  }
  ## ----------------------------------------------------------
  ## compile facs plot
  ## ----------------------------------------------------------

  if (length(facs.plot) != 0) {

    plot <- wrap_plots(facs.plot)

  } else
    plot <- NULL

  ## ----------------------------------------------------------
  ## compile percentage results
  ## ----------------------------------------------------------

  compiled.results <-
    data.frame(variable = character(), percent = numeric())

  for (r in 1:nrow(master.table)) {
    if (master.table[r, "appliedGate"] == "Root") {
      if (!is.na(master.table[r, "qoi"])) {
        df <- data.list[[toupper(paste0(master.table[r, "appliedGate"], "_",
                                        master.table[r, "qoi"]))]]
      } else {
        df <- data.list[[toupper(paste0(master.table[r, "appliedGate"]))]]
      }
      total.cells <- nrow(df)
    } else if (master.table[r, "appliedGate"] != "Root") {
      if (master.table[r, "markerDistribution"] %in% "*,*") {
        lists <-
          data.list[grepl(toupper(paste0(master.table[r, "createdGate"])),
                          names(data.list))]
        df <- plyr::ldply(lists, data.frame)
      } else if (!is.na(master.table[r, "qoi"])) {
        df <- data.list[[toupper(paste0(master.table[r, "createdGate"], "_",
                                        master.table[r, "qoi"]))]]
        if (is.null(df)) {
          df <- data.list[[toupper(paste0(master.table[r, "createdGate"]))]]
        }
      } else {
        df <- data.list[[toupper(paste0(master.table[r, "createdGate"]))]]
      }

      percent.total <-
        nrow(df[df[["phenotypes"]] %in% master.table[r, "phenotype"],]) * 100 /
        total.cells

      compiled.results <-
        as.data.frame(rbind(as.matrix(compiled.results),
                            c(
                              paste0("percent.gated_", master.table[r, "phenotype"]),
                              percent.total
                            )))

      # percent.gate <- nrow(df[df[["phenotypes"]] %in% master.table[r,"phenotype"],])*100 /
      #   nrow(df)
      #
      # compiled.results <- as.data.frame(rbind(as.matrix(compiled.results),
      #                                         c(paste0("percent.gate_",master.table[r,"phenotype"]), percent.gate)))
    }
  }
  ### tidy results

  compiled.results[["percent"]] <-
    as.numeric(as.character(compiled.results[["percent"]]))
  ## ----------------------------------------------------------
  ## merge list to get exact phenotypes
  ## ----------------------------------------------------------

  cordinate.data <- data.list[["ROOT"]]

  df <-
    cordinate.data[, colnames(cordinate.data) %in% "phenotypes", drop = FALSE]


  for (i in 1:length(data.list)) {
    df_uni <-
      data.list[[i]][, colnames(data.list[[i]]) %in% "phenotypes", drop = FALSE]

    df <- merge(df, df_uni, by = 'row.names', all = TRUE) %>%
      column_to_rownames("Row.names")
  }

  df <- as.data.frame(apply(df, 2, function(x)
    gsub('\\s+', '_', x)))

  ## collapse columns
  df[["phenotypes"]] <-
    apply(df, 1, function(x)
      paste(x[!is.na(x)], collapse = "|"))

  ## Select the phenotypes
  df[["phenotypes"]] <-
    sapply(strsplit(df[["phenotypes"]], split = "|", fixed = TRUE), tail, 1L)

  ## tidy phenotypes
  df[["phenotypes"]] <-
    ifelse(df[["phenotypes"]] %in% "character(0)", NA_character_,
           df[["phenotypes"]])



  ## finalize cordinate data
  try(cordinate.data[["phenotypes"]] <-
        unlist(df[["phenotypes"]][match(rownames(cordinate.data),
                                        rownames(df), nomatch = 0)]), silent = TRUE)
  ## ----------------------------------------------------------
  ## add total cell  counts
  ## ----------------------------------------------------------
  for (j in unique(cordinate.data$phenotypes)) {
    percent.total <-
      nrow(cordinate.data[cordinate.data[["phenotypes"]] %in% j,]) * 100 /
      total.cells

    compiled.results <-
      as.data.frame(rbind(as.matrix(compiled.results),
                          c(
                            paste0("percent.total_", j), percent.total
                          )))

  }

  return(
    list(
      cordinate.data = cordinate.data,
      plot = plot,
      summary.results = compiled.results
    )
  )


}
