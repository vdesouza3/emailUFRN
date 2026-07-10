# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^dez.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes

# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_dez.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^nov.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes

# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_nov.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^out.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes

# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_out.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())


###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^set.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes

# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_set.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^ago.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes

# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_ago.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^jul.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes

# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_jul.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^jun.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes


# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_jun.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^mai.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes


# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_mai.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^abr.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes


# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_abr.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^mar.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes


# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_mar.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^fev.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes

# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_fev.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())


###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Carregar pacotes
library(ncdf4)
library(raster)
library(metR)
library(ggplot2)
library(terra)
library(dplyr)
library(sf)

# Caminho para os arquivos
caminho <- "/home/laico/Área de Trabalho/VRSS/dados selecionados"
arquivos <- list.files(caminho, pattern = "^jan.*\\.nc$", full.names = TRUE)

# Diretório para salvar imagens
dir_imagens <- "/home/laico/Área de Trabalho/VRSS/imagens_vetoriais/"
dir.create(dir_imagens, showWarnings = FALSE, recursive = TRUE)

# Carregar shapefiles
shapefile_paises <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/GEOFT_PAIS.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)
shapefile_estados <- st_read("/home/laico/Área de Trabalho/VRSS/dados selecionados/shap/BR_UF_2024.shp", quiet = TRUE) |> 
  st_transform(crs = 4326)

# Alturas e nomes desejados
alturas <- c(120, 150, 180)
nomes_u <- c("cemdoisu", "cemcincou", "cemoitou")
nomes_v <- c("cemdoisv", "cemcincov", "cemoitov")

# Definir anos e regiões de interesse
anos_desejados <- c(2007)  # <- MODIFIQUE AQUI os anos desejados
regioes <- list(                      # <- MODIFIQUE AQUI as regiões (lon_min, lon_max, lat_min, lat_max)
  c(-54, -27, -20, -11)
)
nomes_regioes <- c("Goiás")  # <- MODIFIQUE AQUI os nomes das regiões, na mesma ordem de regioes

# Função para interpolar u e v
interpolar_uv <- function(u10, v10, u100, v100, z) {
  u_z <- u10 + (log(z/10) / log(100/10)) * (u100 - u10)
  v_z <- v10 + (log(z/10) / log(100/10)) * (v100 - v10)
  list(u = u_z, v = v_z)
}

message("Script iniciado em: ", Sys.time())

for (arquivo in arquivos) {
  nc <- nc_open(arquivo)
  anos_todos <- 1961:2024
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  for (i in which(anos_todos %in% anos_desejados)) {
    # carrega componentes do ano escolhido
    u10   <- ncvar_get(nc, "10u",  start = c(1,1,i), count = c(-1,-1,1))
    v10   <- ncvar_get(nc, "10v",  start = c(1,1,i), count = c(-1,-1,1))
    u100  <- ncvar_get(nc, "100u", start = c(1,1,i), count = c(-1,-1,1))
    v100  <- ncvar_get(nc, "100v", start = c(1,1,i), count = c(-1,-1,1))
    
    for (j in seq_along(alturas)) {
      h       <- alturas[j]
      interp  <- interpolar_uv(u10, v10, u100, v100, h)
      vel     <- sqrt(interp$u^2 + interp$v^2)
      
      # con‑verte para data‑frame
      df <- expand.grid(lon = lon, lat = lat) |>
        mutate(
          vel = as.vector(vel),
          u   = as.vector(interp$u),
          v   = as.vector(interp$v),
          ix  = round(lon, 1),
          iy  = round(lat, 1)
        ) |> filter(!is.na(vel))
      
      # ---------- LOOP DE REGIÕES ----------
      for (k in seq_along(regioes)) {
        regiao      <- regioes[[k]]
        nome_regiao <- nomes_regioes[k]
        lon_min <- regiao[1]; lon_max <- regiao[2]
        lat_min <- regiao[3]; lat_max <- regiao[4]
        
        df_recorte <- df |> 
          filter(lon >= lon_min, lon <= lon_max,
                 lat >= lat_min, lat <= lat_max)
        
        # amostra vetores (menos flechas)
        df_setas <- df_recorte |>
          distinct(ix, iy, .keep_all = TRUE) |>
          sample_frac(0.14) |>
          mutate(
            xend = lon + u * 0.09,
            yend = lat + v * 0.09
          ) |>
          filter(xend >= lon_min, xend <= lon_max,
                 yend >= lat_min, yend <= lat_max)
        
        p <- ggplot() +
          geom_raster(data = df_recorte,
                      aes(lon, lat, fill = vel)) +
          scale_fill_viridis_c(
            name   = NULL,
            option = "D",
            limits = c(0, 35),
            breaks = seq(0, 35, 5),
            guide  = guide_colorbar(barwidth = 0.5, barheight = 10,
                                    frame.col = "black", ticks.col = "black")
          ) +
          geom_segment(data = df_setas,
                       aes(x = lon, y = lat, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.08, "cm")),
                       colour = "black", linewidth = 0.3) +
          geom_sf(data = shapefile_paises,  colour = "black", fill = NA, linewidth = 0.3) +
          geom_sf(data = shapefile_estados, colour = "black", fill = NA, linewidth = 0.2) +
          coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
          labs(title = sprintf("Velocidade do vento a %dm (%s de %d) - %s",
                               h, toupper(substr(basename(arquivo),1,4)), anos_todos[i], nome_regiao),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 9) +
          theme(plot.background = element_rect(fill = "white", colour = NA),
                legend.position   = "right",
                legend.justification = c("left", "center"),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
        
        img_nome <- sprintf("vento_%dm_%d_%s_jan.png", h, anos_todos[i], nome_regiao)
        ggsave(filename = file.path(dir_imagens, img_nome), plot = p, width = 8, height = 6)
        cat(sprintf("✅ Imagem salva: %s\n", img_nome))
      }
    }
  }
  nc_close(nc)
  message("✔️ Arquivo concluído:", basename(arquivo))
}

message("\n✅ Processamento completo em: ", Sys.time())

