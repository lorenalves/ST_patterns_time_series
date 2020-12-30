## Exploring spatiotemporal patterns in Cerrado Biome
library(sits)
library(kohonen)
library(dplyr)
library(ggplot2)
library(factoextra)
library(magrittr)


#Call functions
source("./functions/point_to_shape.R")
source("./functions/functions_ST_patterns.R")

# -------------------------------------- sits ---------------------------------------

set.seed(1000)
som_cluster_crop_pasture.lst <-
  sits::sits_som_map(
    samples_crop_pasture.tb,
    grid_xdim = 14,
    grid_ydim = 14,
    alpha = c(0.5,0.1),
    rlen = 100,
    distance = "euclidean"
  )

#only renaming the variable
som_cluster <- som_cluster_crop_pasture.lst

#(Change whatmap by band in sits package)
plot(som_cluster, whatmap = 1)

#Confusion Matrix, purity
evaluate_cluster <- sits_som_evaluate_cluster(som_cluster)
plot(evaluate_cluster, "Cluster Purity")


#From here it is necessarsy call this script
# ---------    Create a new way to plot in SOM map   --------------------
sits_plot_som_map_V2 (som_cluster, type= "by_year", class = "Cropland")

# add neuron label an id_neuron information (it must be internal fuction)
samples_output_som.tb <- .summary_som_output(som_cluster)

# ----------------------------------- Analysis by class -----------------------
class <- "Cropland"

How_many_neurons <- length(which(som_cluster$som_properties$neuron_label == class)) 

#Apply hierarchical clustering on the weight vectors by class


crop.sub <-
  apply_hc_cluster(som_cluster, class = class, linkage = "average")
plot_dendrogram (crop.sub, k =9)

# ------------------------------== Explore clusters --------------------------------
#source("./functions/ST_analysis_V2.R")

n_cluster <-9
get_sub_samples <- info_subgroup(som_cluster, crop.sub, n_cluster)

#change axis labels
plot_info_samples(get_sub_samples)
#plot_info_samples_plotly(get_sub_samples)

#Create a generic plot
plot_subgroups(som_cluster, crop.sub, n_cluster,type_plot = "som")
plot_subgroups(som_cluster, crop.sub, n_cluster,type_plot = "ts", ts_value = "median")
plot_subgroups(som_cluster, crop.sub, n_cluster,type_plot = "ts", ts_value = "together")
plot_subgroups(som_cluster, crop.sub, n_cluster,type_plot = "single_group", group = 3)


plot_sample_by_year (som_cluster, get_sub_samples, subgroups_properties, class = NULL)

#plot paper
plot_subgroups_normal2(som_cluster, crop.sub, n_cluster,type_plot = "ts", ts_value = "together")

# Fix this function
# select_neurons_to_plot(
#   som_cluster,
#   hc = crop.sub,
#   number_of_cluster = n_cluster,
#   neurons =  c(12,34,101)
#   
# )

# -------------------------------- Explore samples ---------------------------------

library(readr)
#Samples west of Bahia
samples_OE_BA_2 <- read_csv("../Tese/samples_OE_BA_2.csv")

#get_sub_samples has the iformation of original time series
plot_unique_sample_by_year (get_sub_samples, samples_OE_BA_2)


# -------------------------------------------------------------------------------------

# Check the number of samples in neuron X
table(dplyr::filter(get_sub_samples, id_neuron %in% c(4,1))$label)

# Check years within a neuron
dplyr::filter(get_sub_samples, id_neuron %in% c(4))
table(get_info_neurons$label, get_info_neurons$year)


#sits plot (original time series)

plot(sits_select_bands(dplyr::filter(get_sub_samples, cluster == 3),ndvi))
plot(sits_select_bands(dplyr::filter(get_sub_samples, id_neuron == 1),ndvi))

# -------------------------------------------------------------------------
library(summarytools) #check this package? Is it necessary?
library(ggpubr)

groups_by_year <- dplyr::select(get_sub_samples, cluster, year)

summarytools::freq(groups_by_year)
contigency_table <- table(year = groups_by_year$year, cluster = groups_by_year$cluster)

ggpubr::ggballoonplot(as.data.frame(contigency_table), fill = "value") +
  ggplot2::scale_fill_viridis_c(option = "C")


# -------------------------------- relabel ----------------------

get_crop <- dplyr::filter(get_sub_samples, label == "Cropland")


Crop_sample <- get_crop %>% mutate(new_label = paste("cluster",cluster, sep="_"))

Crop_sample$label <- Crop_sample$new_label

## -----------------------------------------------------------------------------
## ------------------------------  Validation  ---------------------------------
## -----------------------------------------------------------------------------


# -------------------------------  Input data ------------------------------

input_data_conf_rfor.tb <-
  sits::sits_kfold_validate(
    Crop_sample,
    folds = 5,
    multicores = 1,
    ml_method = sits_rfor()
  )

print("== Confusion Matrix = RFOR =======================")
input_conf_rfor.mx <- sits_conf_matrix(input_data_conf_rfor.tb)

# -------------------------------  Output data ------------------------------
output_data_conf_rfor.tb <-
  sits::sits_kfold_validate(
    output_data_set.tb,
    folds = 5,
    multicores = 1,
    ml_method = sits_rfor()
  )

print("== Confusion Matrix = RFOR =======================")
output_conf_rfor.mx <- sits_conf_matrix(output_data_conf_rfor.tb)

#------------------------------------
