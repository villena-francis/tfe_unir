################################################################################
# SECTION 0: PREPROCESSING - CONVERT TSV FILES TO COMBINED CSV
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Function to extract information from sample ID
extract_sample_info <- function(sid) {
  # Expected format: normal_40x_30000 or tumor_40x_15000
  parts <- strsplit(sid, "_")[[1]]
  
  sample_type <- parts[1]  # normal or tumor
  sequencing_depth <- as.numeric(gsub("x", "", parts[2]))  # 40x -> 40
  read_length <- as.numeric(parts[3])  # 15000, 30000, etc.
  
  return(c(sample_type, sequencing_depth, read_length))
}

# Function to parse CPU Efficiency (e.g., "46.72% of 54-02:48:04 core-walltime")
parse_cpu_efficiency <- function(efficiency_str) {
  # Extract percentage
  percentage <- as.numeric(gsub("%.*", "", efficiency_str))
  return(percentage)
}

# Function to parse Memory Utilized (e.g., "139.82 GB")
parse_memory <- function(memory_str) {
  # Extract numeric value
  memory_gb <- as.numeric(gsub(" GB", "", memory_str))
  return(memory_gb)
}

# Function to process a single TSV file
process_visor_tsv <- function(filepath) {
  # Read TSV file
  data <- read.delim(filepath, stringsAsFactors = FALSE)
  
  # Extract sample information
  sample_info <- t(sapply(data$sid, extract_sample_info))
  
  # Create processed dataframe
  processed_data <- data.frame(
    Sample_type = sample_info[, 1],
    Sequencing_depth = as.numeric(sample_info[, 2]),
    Read_length = as.numeric(sample_info[, 3]),
    CPU_asigned = data$`Cores.per.node`,
    CPU_usage = sapply(data$`CPU.Efficiency`, parse_cpu_efficiency),
    Time = data$`Job.Wall.clock.time`,
    Memory_asigned_GB = 360,
    RAM_usage = sapply(data$`Memory.Utilized`, parse_memory),
    stringsAsFactors = FALSE
  )
  
  return(processed_data)
}

# Process both TSV files
normal_data <- process_visor_tsv("data/cluster_bmk/hpc_data/visor_laser_normal.tsv")
tumor_data <- process_visor_tsv("data/cluster_bmk/hpc_data/visor_laser_tumour.tsv")

# Combine datasets
combined_data <- rbind(normal_data, tumor_data)

# Save to CSV
write.csv(combined_data, "hpc_visor_laser.csv", row.names = FALSE)

################################################################################
# SECTION 1: VISOR ANALYSIS
################################################################################

# Read VISOR performance data
visor_data <- read.csv("hpc_visor_laser.csv")
visor_data$Read_length <- factor(visor_data$Read_length)

# Time conversion function: converts HH:MM:SS to hours
convert_time_to_hours <- function(time_str) {
  # Separate days if they exist
  parts <- strsplit(time_str, "-")[[1]]
  
  if (length(parts) == 2) {
    days <- as.numeric(parts[1])
    time <- parts[2]
  } else {
    days <- 0
    time <- parts[1]
  }
  
  # Convert HH:MM:SS format
  time_components <- as.numeric(strsplit(time, ":")[[1]])
  hours <- time_components[1]
  minutes <- time_components[2]
  seconds <- time_components[3]
  
  # Convert everything to hours
  total_hours <- days * 24 + 
    hours + 
    minutes/60 + 
    seconds/3600
  
  return(round(total_hours, 2))
}

# Process time data and calculate CPU cores
visor_data$Time <- sapply(visor_data$Time, convert_time_to_hours)
visor_data <- visor_data %>%
  mutate(CPU_cores = CPU_asigned * (CPU_usage / 100),
         Read_length_Mb = as.numeric(as.character(Read_length)) / 1000)

# Convert Read_length_Mb to factor for plotting
visor_data$Read_length_Mb <- factor(visor_data$Read_length_Mb)

# Define color schemes for visualization UNIR color based #0098cd
fill_colors <- c(
  CPU_cores = "#0098cd",      
  Time = "#80cbE6",           
  RAM_usage = "#006C8F"       
)

point_colors <- c(
  CPU_cores = "#006D92",      
  Time = "#4DA8C7",           
  RAM_usage = "#004D67"      
)

# Calculate statistics for plotting
avg_CPU_asigned <- mean(visor_data$CPU_asigned)
global_y_max <- max(visor_data$CPU_cores)

# Create CPU usage plot
# Create CPU usage plot
visor_cpu <- ggplot(visor_data, aes(x = Read_length_Mb, y = CPU_cores)) +
  geom_boxplot(
    fill = fill_colors["CPU_cores"],
    color = point_colors["CPU_cores"],
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = point_colors["CPU_cores"],
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  labs(
    x = "Longitud de lectura (Kb)",
    y = "Núcleos CPU"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

# Create memory usage plot
visor_mem <- ggplot(visor_data, aes(x = Read_length_Mb, y = RAM_usage)) +
  geom_boxplot(
    fill = fill_colors["RAM_usage"],
    color = point_colors["RAM_usage"],
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = point_colors["RAM_usage"],
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  labs(
    x = "Longitud de lecturas (Kb)",
    y = "Uso de RAM (GB)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

# Create execution time plot
global_y_max <- max(visor_data$Time)

visor_times <- ggplot(visor_data, aes(x = Read_length_Mb, y = Time)) +
  geom_boxplot(
    fill = fill_colors["Time"],
    color = point_colors["Time"],
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = point_colors["Time"],
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  labs(
    x = "Longitud de lecturas (Kb)",
    y = "Tiempo (horas)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

# Combine VISOR plots
visor_combined_plot <- (visor_cpu + plot_spacer() + visor_mem + plot_spacer() + visor_times) +
  plot_layout(
    widths = c(1, 0.1, 1, 0.1, 1),
    guides = "collect"
  ) &
  theme(aspect.ratio = 1)

print(visor_combined_plot)


################################################################################
# SECTION 2: SV CALLERS ANALYSIS (Adapted for TSV files)
################################################################################

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)

# Function to read and process a TSV file
read_caller_tsv <- function(file_path, caller_name, cpu_assigned, memory_assigned) {
  
  # Read the TSV file
  data <- read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Process the data
  processed_data <- data %>%
    mutate(
      Caller = caller_name,
      # Extract read length from sid (e.g., "v1_40x_15000_B" -> 15 Kb)
      Read_length = as.numeric(sub(".*40x_(\\d+)_.*", "\\1", sid)) / 1000,
      CPU_asigned = cpu_assigned,
      # Parse CPU efficiency percentage (e.g., "70.95%" -> 70.95)
      CPU_usage = as.numeric(sub("%.*", "", `CPU.Efficiency`)),
      # Parse time in HH:MM:SS format and convert to minutes
      Time = sapply(`Job.Wall.clock.time`, function(t) {
        parts <- as.numeric(unlist(strsplit(as.character(t), ":")))
        return((parts[1] * 3600 + parts[2] * 60 + parts[3]) / 60)
      }),
      Memory_asigned = memory_assigned,
      # Parse memory (e.g., "61.86 GB" -> 61.86)
      RAM_usage = as.numeric(sub(" GB", "", `Memory.Utilized`))
    ) %>%
    select(Caller, Read_length, CPU_asigned, CPU_usage, Time, Memory_asigned, RAM_usage)
  
  return(processed_data)
}

# Read all TSV files
savana_data <- read_caller_tsv(
  "data/cluster_bmk/hpc_data/savana.tsv", 
  "SAVANA", 
  cpu_assigned = 24, 
  memory_assigned = 400
)

severus_data <- read_caller_tsv(
  "data/cluster_bmk/hpc_data/severus.tsv", 
  "Severus", 
  cpu_assigned = 24, 
  memory_assigned = 40
)

sniffles_data <- read_caller_tsv(
  "data/cluster_bmk/hpc_data/sniffles.tsv", 
  "Sniffles2", 
  cpu_assigned = 24, 
  memory_assigned = 40
)

svision_data <- read_caller_tsv(
  "data/cluster_bmk/hpc_data/svision.tsv", 
  "Svision-pro", 
  cpu_assigned = 52, 
  memory_assigned = 40
)

# Combine all data
hpc_data <- bind_rows(savana_data, severus_data, sniffles_data, svision_data)

# Convert Read_length to factor
hpc_data$Read_length <- factor(hpc_data$Read_length)

# Calculate CPU cores
hpc_data <- hpc_data %>%
  mutate(CPU_cores = CPU_asigned * (CPU_usage / 100))

# Transform data to long format
hpc_data_long <- hpc_data %>%
  pivot_longer(
    cols = c(CPU_cores, RAM_usage, Time),
    names_to = "Metric", 
    values_to = "Value"
  )

# Calculate CPU statistics
cpu_assigned_avg <- hpc_data %>%
  group_by(Caller) %>%
  summarise(avg_CPU_asigned = mean(CPU_asigned, na.rm = TRUE))

# Join average CPU data
hpc_data_long <- hpc_data_long %>%
  left_join(cpu_assigned_avg, by = "Caller")

# Calculate global maximum for CPU cores
global_y_max <- hpc_data_long %>%
  filter(Metric == "CPU_cores") %>%
  summarise(global_y_max = max(Value, na.rm = TRUE)) %>%
  pull(global_y_max)

# Define color schemes for visualization UNIR color based #0098cd
fill_colors <- c(
  CPU_cores = "#0098cd",
  Time = "#80cbE6",
  RAM_usage = "#006C8F"
)

point_colors <- c(
  CPU_cores = "#006D92",
  Time = "#4DA8C7",
  RAM_usage = "#004D67"
)

# Function to darken colors
darken <- function(color, factor = 0.2) {
  col_rgb <- col2rgb(color) / 255
  darkened <- rgb(col_rgb[1] * (1 - factor), 
                  col_rgb[2] * (1 - factor), 
                  col_rgb[3] * (1 - factor))
  return(darkened)
}

# Crear etiquetas en español con unidades
metric_labels <- c(
  CPU_cores = "Núcleos CPU",
  RAM_usage = "Uso de RAM (GB)",
  Time = "Tiempo (minutos)"
)

# Create SV callers performance plot
hpc_plot <- ggplot(hpc_data_long, aes(x = Read_length, y = Value)) +
  geom_boxplot(
    aes(fill = Metric, color = Metric),
    alpha = 0.5
  ) +
  geom_jitter(
    aes(color = Metric),
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  geom_hline(
    data = subset(hpc_data_long, Metric == "CPU_cores"), 
    aes(yintercept = avg_CPU_asigned),
    linetype = "dashed",
    color = "darkgray"
  ) +
  geom_text(
    data = subset(hpc_data_long, Metric == "CPU_cores"),
    aes(
      x = Inf,
      y = avg_CPU_asigned,
      label = paste("asignado =", round(avg_CPU_asigned, 2)),
      vjust = ifelse(avg_CPU_asigned > (global_y_max * 0.9), 1.5, -0.5)
    ),
    hjust = 1.1,
    color = "darkgray",
    size = 3
  ) +
  facet_grid(rows = vars(Metric), cols = vars(Caller), scales = "free_y",
             labeller = labeller(Metric = metric_labels)) +
  scale_x_discrete(breaks = c("15", "30", "50", "100")) +
  labs(x = "Longitud de lecturas (Kb)", y = "") +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(
    values = c(point_colors, sapply(fill_colors, function(color) darken(color, 0.2)))
  )

print(hpc_plot)

# Optional: Save the combined data to verify the transformation
write.csv(hpc_data, "hpc_data_from_tsv.csv", row.names = FALSE)


################################################################################
# SECTION 3: CALLER ACCURACY ANALYSIS
################################################################################

# Read accuracy metrics data
calls_data <- read.csv("data/cluster_bmk/calls_data/callings.csv")

# Rename columns to match expected format
calls_data <- calls_data %>%
  rename(
    Caller = caller,
    Read_length = read_length
  )

# Convert read_length from bp to Kb and make it a factor
calls_data <- calls_data %>%
  mutate(Read_length = Read_length / 1000)

calls_data$Read_length <- factor(calls_data$Read_length)

# Transform to long format
calls_data_long <- calls_data %>%
  pivot_longer(
    cols = c(Precision, Recall, F1_score),
    names_to = "Metric",
    values_to = "Value"
  )

# Order metric levels
calls_data_long$Metric <- factor(
  calls_data_long$Metric, 
  levels = c("Precision", "Recall", "F1_score")
)

# Define colors for accuracy metrics using UNIR palette
fill_colors_accuracy <- c(
  Precision = "#0098cd",      
  Recall = "#80cbE6",         
  F1_score = "#006C8F"        
)

point_colors_accuracy <- c(
  Precision = "#006D92",      
  Recall = "#4DA8C7",         
  F1_score = "#004D67"        
)

# Crear etiquetas en español
metric_labels_accuracy <- c(
  Precision = "Precisión",
  Recall = "Sensibilidad",
  F1_score = "F1-score"
)

# Create accuracy metrics plot
calls_plot <- ggplot(calls_data_long, aes(x = Read_length, y = Value)) +
  geom_boxplot(
    aes(fill = Metric, color = Metric),
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = Metric),
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  facet_grid(rows = vars(Metric), cols = vars(Caller), scales = "free_y",
             labeller = labeller(Metric = metric_labels_accuracy)) +
  scale_x_discrete(breaks = c("15", "30", "50", "100")) +
  labs(x = "Longitud de lecturas (Kb)", y = "") +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = fill_colors_accuracy) +
  scale_color_manual(values = point_colors_accuracy)

print(calls_plot)
