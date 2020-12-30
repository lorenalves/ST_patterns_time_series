#Join the tibbles in output of SOM
# Suggest to Gilberto include it in sits (think how to present it)
.summary_som_output <- function (som_cluster){
  
  # Get id_sample, id_neuron and neuron_label from som_cluster
  id_neuron_sample <- unique(dplyr::select(som_cluster$statistics_samples$samples_t, id_sample, id_neuron, neuron_label))
  
  #som_output + id_neuron
  samples_complete <- som_cluster$samples_output.tb %>% inner_join(id_neuron_sample)
  
  return (samples_complete)
  
}

# ------------------------------------------------------------------------------
#
#' @title Apply hierarchical clustering
#' @name apply_hc_cluster
#' @keywords Hiearchical clustering
#' @author Lorena Santos, \email{lorena.santos@@inpe.br}
#'
#' @description        Apply hierarchical clustering over the weight vectors
#' @param som_cluster  Object returned by sits_som_cluster
#' @param class        Class label
#' @param linkage        Linkage criteria (ward.D2, single, centroids, average)
#' @return A object from hclust (stats) that describe the tree produced by HC.
#' ------------------------------------------------------------------------------

apply_hc_cluster <- function (som_cluster, class = NULL,  linkage = "ward.D2")
{
  #filter the samples by neuron_label from class_neurons
  samples_information <-
    dplyr::filter(som_cluster$statistics_samples$samples_t,
                   neuron_label == class)
  
  
  
  #get the neuron's id of where these samples were allocated
  neurons_class <- unique(samples_information$id_neuron)
  
  #get the  weigth's vector of neurons
  codebooks <- som_cluster$som_properties$codes
  
  #*********  It must be dynamic *******
  #get the ndvi weight
  # weight_ndvi <- som_cluster$som_properties$codes$ndvi
  # weight_evi <- som_cluster$som_properties$codes$evi
  # weight_mir <- som_cluster$som_properties$codes$mir
  # weight_nir <- som_cluster$som_properties$codes$nir
  # 
  # all_indices <- cbind(weight_ndvi, weight_evi, weight_nir, weight_mir)
  
  all_indices <- do.call(cbind, codebooks)
  
  
  # Only neurons of class X
  weight.ts <- all_indices[neurons_class, ]
  
  
  #(future - Create a parameter to choose the distance )
  distance_atrributes <- proxy::dist(weight.ts, distance = "euclidean")
  hc <- stats::hclust(distance_atrributes, linkage)
  
  return (hc)
  
} 

# ------------------------------------------------------------------------------
#
#' @title Plot dendrogram
#' @name plot_dendrogram
#' @keywords dendrogram
#' @author Lorena Santos, \email{lorena.santos@@inpe.br}
#'
#' @description       Create a dendrogram based on package factoextra (enhanced dendrogram)
#' @param hc          Object returned hierarchical clustering
#' @param k           Set number of clusters
#' @param h           Cut the dendrogram in height h

#' ------------------------------------------------------------------------------
plot_dendrogram <- function (hc, k = NULL, h = NULL)
{
  plot_dendro <- factoextra::fviz_dend(
    hc,
    k = k, # Cut in k groups
    h = h, # Cut in height h
    cex = 0.6, # label size
    rect = TRUE, # Add rectangle around groups
    rect_fill = TRUE, 
    rect_border = "jco",#color
    color_labels_by_k = TRUE # color labels by groups
    
  )
  return (plot_dendro)
  
}

# ------------------------------------------------------------------------------

info_subgroup <- function (som_cluster, hc, number_of_cluster)
{
  
  info_cluster <- .get_info_subgroups (som_cluster, hc, number_of_cluster)
  samples_output_som.tb <- .summary_som_output(som_cluster)
  
  #get samples assigned to neurons labeled as class X
  get_samples <- dplyr::filter(samples_output_som.tb, id_neuron %in% info_cluster$cluster_neuron$id_neuron)
 
  get_sub_samples <-
    get_samples %>% dplyr::inner_join(info_cluster$cluster_neuron) %>% dplyr::mutate(year = lubridate::year(start_date))
  
  return (get_sub_samples)
  
}

plot_subgroups <- function (som_cluster, hc, number_of_cluster, clean_color_neurons = TRUE, ts_value = "together", type_plot ="ts", group = NULL)
{
 
  info_cluster <- .get_info_subgroups (som_cluster, hc, number_of_cluster, ts_value)
  
  #create a unique object to be a parameter of plots
  
  if (type_plot == "ts"){
    .ts_subgroups_plot(info_cluster$data_frame_info_cluster,info_cluster$label_class, ts_value)
  }else if (type_plot == "som")
  {
    .som_subgroups_plot(som_cluster, number_of_cluster,info_cluster$cluster_neuron,info_cluster$label_class, clean_color_neurons = FALSE)
  } else if (type_plot == "single_group")
  {
     info_cluster <- .get_info_subgroups (som_cluster, hc, number_of_cluster, ts_value = "single")
    .ts_one_subcluster(info_cluster$data_frame_info_cluster, group)
  }
  #return (invisible(cluster_neuron))
  
}

# Fixed point and change in time
plot_unique_sample_by_year <- function(get_sub_samples, samples_same_point)
{
  sample_same_SPL <- dplyr::filter(get_sub_samples, id_sample %in% samples_same_point$id_sample)
  
  #prepare  data frame for plotting
  plot.df <- data.frame()
  sample_same_SPL$pattern <- paste(sample_same_SPL$year, " - Neuron ", sample_same_SPL$id_neuron, sep = "")
  data.tb <- sample_same_SPL
  
  lat <- unique(data.tb$latitude)
  long <- unique(data.tb$longitude)
  
  
  # put the time series in the data frame
  purrr::pmap(list(data.tb$pattern, data.tb$time_series),
              function(pattern, ts) {
                lb <- as.character(pattern)
                # extract the time series
                df <- data.frame(Time = ts$Index, ts[-1], Pattern = lb)
                plot.df <<- rbind(plot.df, df)
              })
  
  plot.df <- reshape2::melt(plot.df, id.vars = c("Time", "Pattern"))
  
  # Plot temporal patterns
  gp <-  ggplot2::ggplot(plot.df, ggplot2::aes_string(x = "Time",
                                                      y = "value",
                                                      colour = "variable") ) +
    ggplot2::geom_line(size = 1) +
    ggplot2::facet_wrap(~Pattern, scales = "free") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_date(labels = scales::date_format("%b")) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Bands")) +
    ggplot2::ggtitle(paste("MODIS time series samples: point (",lat, ",",long,")", sep = " "))+
    ggplot2::ylab("Value")+
    ggplot2::coord_cartesian(ylim = c(0, 1))
   
  plotly::ggplotly() 
  #graphics::plot(gp)
}

#BarPlot (plot the number of samples by cluster)
plot_info_samples <- function (get_sub_samples)
{
  #get id_sample, label, cluster and id_neuron
  samples.temp <- dplyr::select(get_sub_samples,id_sample, label, cluster, id_neuron)
  
  #get cluster and id_neuron
  cluster_neuron_id <- unique(dplyr::select(get_sub_samples, cluster, id_neuron))
  
  #How many samples (by class) there are in each cluster?
  info_samples <- samples.temp %>% dplyr::group_by(cluster,label) %>% dplyr::summarise(n_samples =  n())
  
  #How many neurons there are in each cluster?
  info_neurons <- cluster_neuron_id %>% dplyr::group_by(cluster) %>% dplyr::summarise(n_neurons =  n())
  
  # Just join the information above
  info_cluster <- tibble::as_tibble(info_samples %>% dplyr::inner_join(info_neurons))
  #cluster <- as.character(info_cluster$cluster)
  
  #samples
  p <- ggplot2::ggplot(data = info_cluster,
                       aes(
                         x = as.factor(info_cluster$cluster),
                         y = info_cluster$n_samples,
                         fill =  info_cluster$label
                       )) +
    ggplot2::geom_bar(position = 'dodge', stat = 'identity') +
    ggplot2::geom_text(
      aes(label = info_cluster$n_samples),
      position = position_dodge(width = 0.9),
      vjust = -0.25
    )+  
    ggplot2::labs(x = "Cluster", y = "Number of samples", colour = "Class") + 
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Number of samples by cluster") +
    ggplot2::scale_fill_discrete(name = "Class")
    #ggplot2::coord_cartesian(ylim = c(0, 1))
  plotly::ggplotly() 
  #graphics::plot(p)
  
}

# -------------------------------- Internal functions -----------------------------------------

.get_info_subgroups <-function (som_cluster, hc, number_of_cluster, ts_value = "together")
{
  df_j.temp <- data.frame()
  df_sub_k.tb <-  data.frame()
  colorful_subgroups <- data.frame()
  subgroups_properties <- NULL
  
  
  #cut dendogram in k clusters
  cut_hc <- stats::cutree(hc, k = number_of_cluster)
  
  #get id_neurons of class X
  id_neurons <- names(cut_hc)
  
  
  # Must be dynamic ( I need to know the attributes here)
  bands_name <- (names(som_cluster$som_properties$codes))
  
  # discover the size of each band
  all_indices <- do.call(cbind, som_cluster$som_properties$codes)
  points_by_band <- dim(all_indices)[2]/length(bands_name)
  
  # weight_ndvi <- som_cluster$som_properties$codes$ndvi
  # weight_evi <- som_cluster$som_properties$codes$evi
  # weight_mir <- som_cluster$som_properties$codes$mir
  # weight_nir <- som_cluster$som_properties$codes$nir
  
  
  #the follow lines do it: weight_ndvi <- som_cluster$som_properties$codes$ndvi
  initial_layer <- 1L
  varname_weights <- paste("weight_",names(som_cluster$som_properties$codes), sep="")
  
  #Get the weight vectors for each band
  for (bn in 1: length(bands_name)){
    
    #all_indice[,1:23] -- dim [169,92]
    end_layer <- bn*points_by_band
    assign(varname_weights[bn], (all_indices[,initial_layer:end_layer]))
    initial_layer <- bn*points_by_band + 1
  }
  
  
  # now, it can be fix..
  index_time <- c("2015-09-14", "2015-09-30", "2015-10-16", "2015-11-01",
                  "2015-11-17", "2015-12-03", "2015-12-19", "2016-01-01", "2016-01-17",
                  "2016-02-02", "2016-02-18", "2016-03-05", "2016-03-21", "2016-04-06",
                  "2016-04-22", "2016-05-08", "2016-05-24", "2016-06-09", "2016-06-25",
                  "2016-07-11", "2016-07-27" ,"2016-08-12","2016-08-28")
  
  
  
  #k is the number of groups
  for (k in 1:length(unique(cut_hc)))
  {
    #get the id_neurons of each subgroups
    neuron_id_string <- which(cut_hc == k)
    neuron_id_string <- names(neuron_id_string)
    neuron_id_int <- as.integer(substring(neuron_id_string, 2))
    
    df_sub_j.tb <-  data.frame()
    
    #build the data.frame by group
    #j is the number of neurons in cluster k
    for (j in 1:length(unique(neuron_id_int))){
      
      
      # It must be on the fly (similar weight_band)
      ndvi= weight_ndvi[neuron_id_int[j],]
      evi = weight_evi[neuron_id_int[j],]
      nir = weight_nir[neuron_id_int[j],]
      mir = weight_mir[neuron_id_int[j],]
      
      
      # same to whatmap in kohonen package
      #position_band <- which(bands_name == band)
      
      if (ts_value == "together")
      {
        #dataframe (id_neuron,  cluster, Time, band)
        df_j.temp <-
          data.frame(
            id_neuron = neuron_id_int[j],
            cluster = k,
            Time = as.Date(index_time),
            ndvi = ndvi #change the name ndvi to selected_band
            #it must be select by a parameter (think about it!)
            
          )
      }else if (ts_value == "median")
      {
        df_j.temp <-
          data.frame(
            #Neuron = neuron_id_int[j],
            cluster = k,
            Time = as.Date(index_time),
            ndvi = ndvi,
            evi = evi,
            nir = nir,
            mir = mir
          )
      }else if (ts_value == "single")
      {
        df_j.temp <-
          data.frame(
            neuron = neuron_id_int[j],
            cluster = k,
            Time = as.Date(index_time),
            ndvi = ndvi,
            evi = evi,
            nir = nir,
            mir = mir
          )
      }
      
      
      #concat lines (neurons)
      df_sub_j.tb <- rbind(df_sub_j.tb, df_j.temp)
      
      # (cluster, id_neuron)
      colorful_subgroups.temp <- data.frame(cluster = k, id_neuron = neuron_id_int)
      
      
    }
    #what neurons must be painted ? (cluster, id_neuron) Full dataframe
    colorful_subgroups <- rbind(colorful_subgroups,colorful_subgroups.temp)
    
    #dataframe (id_neuron, cluster, time, selected_band)
    df_sub_k.tb <-rbind (df_sub_k.tb, df_sub_j.tb)
  }
  
  #tibble (cluster, id_neuron)
  colors_subgroups <- tibble::as_tibble(colorful_subgroups)
  cluster_neuron <- colors_subgroups
  #what is the class?
  label_class.tb <-
    unique(dplyr::filter(som_cluster$statistics_samples$neuron_t,
                         id_neuron == neuron_id_int[1])$neuron_label)
  
  subgroups.lst <-
    structure(
      list(
        data_frame_info_cluster = df_sub_k.tb,
        cluster_neuron = cluster_neuron,
        label_class = as.character(label_class.tb)
      )
    )
  
  return(subgroups.lst)
  
  
}

#Plot clusters
.ts_subgroups_plot <- function (df_sub_k.tb , label_class,ts_value){
  
  if (ts_value == "together")
  {
    #All plots - one band
    # Melt (time, cluster, id_neuron, variable (name band), value of band)
    df2 <- reshape2::melt(df_sub_k.tb, id.vars = c("Time", "cluster", "id_neuron"))
    df3 <- tibble::as_tibble(df2)
    df3$variable <- as.character(df3$id_neuron)
  }else if (ts_value == "median")
  {
    df2 <- reshape2::melt(df_sub_k.tb, id.vars = c("Time", "cluster"))
    df3 <- tibble::as_tibble(df2)
    df3 <- df3 %>%
      select(Time, cluster, variable, value) %>%
      group_by(Time, cluster, variable) %>%
      summarise(value  = median(value))
  }
  
  # Plot median of the groups
  gp <-  ggplot2::ggplot(df3, ggplot2::aes_string(x = "Time",
                                                  y = "value",
                                                  colour = "variable")) +
    ggplot2::geom_line(size = 0.3) +
    ggplot2::facet_wrap(~cluster, scales = "free") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_x_date(labels = scales::date_format("%b")) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Bands")) +
    ggplot2::ylab("Values")+
    ggplot2::ggtitle(paste("Subgroups of ",label_class, sep = " ")) +
    ggplot2::stat_summary(fun.data ="mean_sdl", geom = "smooth")+
    ggplot2::coord_cartesian(ylim = c(0, 1))
  plotly::ggplotly() 
  #ggplot2::theme_minimal()
  #graphics::plot(gp)
}

#plot som map
.som_subgroups_plot <- function (som_cluster,number_of_cluster, colors_subgroups, label_class, clean_color_neurons = TRUE){
  #Plot sugbroups on SOM GRID
  Dark2 <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
  BrBG <- RColorBrewer::brewer.pal(n = 8, name = "BrBG")
  paired <- RColorBrewer::brewer.pal(n = 8, name = "Paired")
  set_colors <- c(Dark2,BrBG,paired)
  
  #cluster and id_neuron
  bgcol_copy2 = som_cluster$som_properties$paint_map
  
  if(clean_color_neurons)
    bgcol_copy2[] <- "white"
  
  #n_cluster <- length(unique(colors_subgroups$cluster))
  for (i in 1: number_of_cluster){
    get_neurons <- dplyr::filter(colors_subgroups, cluster == i)$id_neuron
    bgcol_copy2[get_neurons] <- set_colors[i]
    
  }
  #get all neurons that belongs to the general class
  #id_neurons_int.sub <- as.integer(substring(id_neurons, 2))
  id_neurons_int.sub <- colors_subgroups$id_neuron
  
  vector_id <- 1:length(bgcol_copy2)
  for(i in 1: length(bgcol_copy2)){
    names(vector_id)[i]<-paste("V",i,sep="")
  }
  vector_id[]<-1
  vector_id[id_neurons_int.sub]<-2
  
  graphics::plot(som_cluster$som_properties,  bgcol = bgcol_copy2 , "codes", whatmap = 1, codeRendering = "lines", main = paste("Subgroups of ",label_class, sep = " "))
  
  #cluster boundaries
  add.cluster.boundaries(som_cluster$som_properties, vector_id)
  
  leg <- unique(paste(label_class ,colors_subgroups$cluster, sep = "_"))
  legend("bottomright", 
         legend = leg, 
         col = set_colors[1:number_of_cluster], 
         pch = c(15), 
         bty = "n", 
         pt.cex = 2, 
         cex = 1, 
         text.col = "black", 
         horiz = F ,
         ncol = 1,
         inset = c(0, 0.025))
}

#Plot all neurons within a cluster
.ts_one_subcluster <- function (df_sub_k.tb, group)
{
  
  df2 <- reshape2::melt(df_sub_k.tb, id.vars = c("Time", "cluster", "neuron"))
  
  df3 <- tibble::as_tibble(df2)
  df3 <- dplyr::filter(df3, cluster == group)
  df3$neuron <- paste("neuron ", df3$neuron, sep = "")
  
  gp <-  ggplot2::ggplot(df3, ggplot2::aes_string(x = "Time",
                                                  y = "value",
                                                  colour = "variable")) +
    ggplot2::geom_line(size = 0.5) +
    ggplot2::facet_wrap(~neuron, scales = "free") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_date(labels = scales::date_format("%b")) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Bands")) +
    ggplot2::ylab("Value")+
    ggplot2::ggtitle(paste("Weight vectors of cluster  ",group, sep = " "))+
    ggplot2::theme_minimal()+
    ggplot2::coord_cartesian(ylim = c(0, 1))
  plotly::ggplotly()
  #graphics::plot(gp)
}

#Fix it
select_neurons_to_plot <- function (som_cluster, hc, number_of_cluster,neurons)
{
  cut_hc <- stats::cutree(hc, k = number_of_cluster)
  id_neurons <- names(cut_hc)
  
  df_sub_k.tb <-  data.frame()
  colorful_subgroups <- data.frame()
  subgroups_properties <- NULL
  
  weight_ndvi <- som_cluster$som_properties$codes$ndvi
  weight_evi <- som_cluster$som_properties$codes$evi
  weight_mir <- som_cluster$som_properties$codes$mir
  weight_nir <- som_cluster$som_properties$codes$nir
  
  index_time <- c("2015-08-29", "2015-09-14", "2015-09-30", "2015-10-16", "2015-11-01",
                  "2015-11-17", "2015-12-03", "2015-12-19", "2016-01-01", "2016-01-17",
                  "2016-02-02", "2016-02-18", "2016-03-05", "2016-03-21", "2016-04-06",
                  "2016-04-22", "2016-05-08", "2016-05-24", "2016-06-09", "2016-06-25",
                  "2016-07-11", "2016-07-27" ,"2016-08-12")
  #k is the number of groups
  for (k in 1:length(unique(cut_hc)))
  {
    #get the id_neurons of each subgroups
    neuron_id_string <- which(cut_hc == k)
    neuron_id_string <- names(neuron_id_string)
    neuron_id_int <- as.integer(substring(neuron_id_string, 2))
    
    df_sub_j.tb <-  data.frame()
    #build the data.frame by group
    for (j in 1:length(unique(neuron_id_int))){
      
      #isso deve ser dinamico no futuro
      ndvi= weight_ndvi[neuron_id_int[j],]
      evi = weight_evi[neuron_id_int[j],]
      nir = weight_nir[neuron_id_int[j],]
      mir = weight_mir[neuron_id_int[j],]
      
      
      df_j.temp <-
        data.frame(
          neuron = neuron_id_int[j],
          cluster = k,
          Time = as.Date(index_time),
          ndvi = ndvi,
          evi = evi,
          nir = nir,
          mir = mir
        )
      
      df_sub_j.tb <- rbind(df_sub_j.tb, df_j.temp)
      colorful_subgroups.temp <- data.frame(cluster = k, id_neuron = neuron_id_int)
      
      
    }
    colorful_subgroups <- rbind(colorful_subgroups,colorful_subgroups.temp)
    df_sub_k.tb <-rbind (df_sub_k.tb, df_sub_j.tb)
  }
  
  label_class <-
    unique(dplyr::filter(som_cluster$statistics_samples$neuron_t,
                         id_neuron == neuron_id_int[1])$neuron_label)
  
  df2 <- reshape2::melt(df_sub_k.tb, id.vars = c("Time", "cluster", "neuron"))
  
  df3 <- as_tibble(df2)
  df3 <- dplyr::filter(df3, neuron %in% neurons)
  df3$neuron <- paste("neuron ", df3$neuron, sep = "")
  
  # 
  # df3 <- df3 %>%
  #   select(Time, cluster, variable, value) %>%
  #   group_by(Time, cluster, variable) %>%
  #   summarise(value  = median(value))
  # # 
  
  # Plot median of the groups
  gp <-  ggplot2::ggplot(df3, ggplot2::aes_string(x = "Time",
                                                  y = "value",
                                                  colour = "variable")) +
    ggplot2::geom_line(size = 1) +
    ggplot2::facet_wrap(~neuron, scales = "free") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_date(labels = scales::date_format("%b")) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Bands")) +
    ggplot2::ylab("Value")+
    #ggplot2::ggtitle(paste("Weight vector of cluster  ",sub_cluster, sep = " "))+
    ggplot2::theme_minimal()+
    ggplot2::coord_cartesian(ylim = c(0, 1))
  
  graphics::plot(gp)
}


#' @title Plot som map
#' @name sits_plot_som_cluster_V2
#' @keywords Suggest sits
#' @author Lorena Santos, \email{lorena.santos@@inpe.br}
#'
#' @description    Plot Som map
#' @param som_cluster      Object returned by sits_som_cluster
#' @param type     type of plot  (codes, mapping, by_year)
#' @param whatmap  band
#' @param class    class label (if the type is by_year)

sits_plot_som_cluster_V2 <- function(som_cluster, type = "codes", whatmap = 1 , class = NULL)
{
  if (type == "mapping") {
    graphics::plot(som_cluster$som_properties,  bgcol = som_cluster$som_properties$paint_map , "mapping", whatmap = whatmap)
  } else if (type == "codes" ){
    graphics::plot(som_cluster$som_properties,  bgcol = som_cluster$som_properties$paint_map , "codes", whatmap = whatmap, codeRendering = "lines")
  }else if (type == "by_year"){
    
    data.tb <- dplyr::select(som_cluster$samples_output.tb,id_sample,latitude,longitude,start_date, end_date,label)
    samples_information <- som_cluster$statistics_samples$samples_t
    it <- unique(max(samples_information$iteration))
    samples_information <- dplyr::filter(samples_information, samples_information$iteration == it)
    samples_st_id <- samples_information %>% dplyr::inner_join(data.tb, by = "id_sample")
    
    
    if (!is.null(class))
    {
      samples_st_id <- dplyr::filter(samples_st_id, samples_st_id$original_label == class)
    }
    
    id_all_year <- samples_st_id %>% dplyr::pull(id_neuron)
    id_all_year <- as.numeric(id_all_year)
    graphics::plot(som_cluster$som_properties,  "mapping", classif = id_all_year , bgcol= som_cluster$som_properties$paint, main = "All years" )
    
    year <- sort(unique(samples_st_id$start_date))
    year <- sort(unique(lubridate::year(samples_st_id$start_date)))
    n_year <- length(year)
    
    for (i in 1:length(year))
    {
      #samples_by_year <- dplyr::filter(samples_st_id, samples_st_id$start_date == year[i])
      
      samples_by_year<- dplyr::filter(samples_st_id, lubridate::year(samples_st_id$start_date) == year[i] )
      text_year <- substr(year[i], 1, 4)
      #get the neuron
      id_samples_year <- samples_by_year %>% dplyr::pull(4)
      id_samples_year <- as.numeric(id_samples_year)
      graphics::plot(som_cluster$som_properties,  "mapping", classif =id_samples_year , bgcol= som_cluster$som_properties$paint, main = text_year)
    }
    
  }
  
  #create a legend (fix it)
  leg <- cbind(som_cluster$som_properties$neuron_label, som_cluster$som_properties$paint_map)
  graphics::legend(
    "bottomright",
    legend = unique(leg[, 1]),
    col = unique(leg[, 2]),
    pch = 15,
    pt.cex = 2,
    cex = 1,
    text.col = "black",
    #horiz = T ,
    inset = c(0.0095, 0.05),
    xpd = TRUE,
    ncol = 1
  )
}

