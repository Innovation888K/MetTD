Load_SBML = function(model_name, option = "human",
                      local_path="C:/Users/梦羽/Desktop/Rspace/MetTD/") {
  library(rsbml)
  options(expressions = 5e5)
  return(list(model = rsbml_read(
    paste(local_path, "SBMLmodels/", model_name, sep = "")
  )))
}

get_header = function(model) {
  com_name = c()
  compartment = c()
  for (i in 1:length(model$model@model@species)) {
    if (length(model$model@model@species[[i]]@name) > 0) {
      com_name[i] = model$model@model@species[[i]]@name
    }
    else{
      com_name[i] = model$model@model@species[[i]]@id
    }
    compartment[i] = model$model@model@compartments[model$model@model@species[[i]]@compartment][[1]]@name
  }
  names = names(model$model@model@species)
  header = data.frame(com_name = com_name,
                      compartment = compartment)
  rownames(header) = names
  model$header = header
  return(model)
}

get_distance_matrix = function(Model) {
  model = Model$model
  header = Model$header
  adjacency_matrix = rep(rep(0, nrow(header)), nrow(header))
  dim(adjacency_matrix) = c(nrow(header), nrow(header))
  adjacency_matrix = data.frame(adjacency_matrix)
  rownames(adjacency_matrix) = rownames(header)
  colnames(adjacency_matrix) = rownames(header)
  pb = txtProgressBar(style = 3)
  index = 1
  for (reaction in model@model@reactions) {
    setTxtProgressBar(pb, index / length(model@model@reactions))
    reac = names(reaction@reactants)
    prod = names(reaction@products)
    modi = names(reaction@modifiers)
    reac_idx = match(reac, rownames(adjacency_matrix))
    prod_idx = match(prod, rownames(adjacency_matrix))
    modi_idx = match(modi, rownames(adjacency_matrix))
    if (length(modi_idx) > 0) {
      if (length(prod_idx) > 0) {
        for (i in 1:length(modi_idx)) {
          for (j in 1:length(prod_idx)) {
            adjacency_matrix[modi_idx[i], prod_idx[j]] = 1
          }
        }
      }
      if (length(reac_idx) > 0) {
        for (i in 1:length(modi_idx)) {
          for (j in 1:length(reac_idx)) {
            adjacency_matrix[modi_idx[i], reac_idx[j]] = 1
          }
        }
      }
    }
    index = index + 1
  }
  close(pb)
  heads = paste(header$com_name, header$compartment, sep = ".")
  freq_matrix = as.data.frame(table(heads))
  remove_symble = freq_matrix$heads[freq_matrix$Freq > 1]
  idx = data.frame(rem = match(remove_symble, heads))
  freq_idx = data.frame(rem = which(heads %in% remove_symble))
  freq_idx = dplyr::anti_join(freq_idx, idx)
  if(length(freq_idx$rem)>0){
    adjacency_matrix = adjacency_matrix[-freq_idx$rem,]
    adjacency_matrix = dplyr::select(adjacency_matrix, colnames(adjacency_matrix)[-freq_idx$rem])
    heads = heads[-freq_idx$rem]
    header = header[-freq_idx$rem,]
  }
  rownames(adjacency_matrix) = heads
  colnames(adjacency_matrix) = heads
  g = igraph::graph_from_adjacency_matrix(as.matrix(adjacency_matrix))
  Model$distance_martix = 2 * data.frame(igraph::distances(g, weights = NA))
  Model$freq_idx = freq_idx
  Model$graph_model=g
  return(Model)
}

get_annotation = function(Model, from = "species") {
  model = Model$model
  annotation = lapply(slot(model@model, from), function(x) {
    if (length(x@cvTerms) > 0) {
      ret = c()
      index = 1
      for (cv in x@cvTerms) {
        for (res in cv@resources) {
          ret[index] = res
          index = index + 1
        }
      }
      return(ret)
    }
  })
  id_styles = c()
  for (i in 1:length(annotation)) {
    for (str in annotation[[i]]) {
      temp = gsub(" ", "", strsplit(str, "/")[[1]])
      if (!temp[4] %in% id_styles) {
        id_styles = c(id_styles, temp[4])
      }
    }
  }
  annotations = rep(rep(NA, length(id_styles)), length(slot(model@model, from)))
  dim(annotations) = c(length(slot(model@model, from)), length(id_styles))
  colnames(annotations) = id_styles
  annotations = data.frame(annotations)
  for (i in 1:length(annotation)) {
    if (!is.null(annotation[[i]])) {
      for (str in annotation[[i]]) {
        temp = gsub(" ", "", strsplit(str, "/")[[1]])
        id_res = gsub(" ", "", strsplit(str, paste(temp[4], "/", sep = ""))[[1]])[2]
        annotations[[temp[4]]][i] = id_res
      }
    }
  }
  if (length(Model$freq_idx$rem) > 0) {
    annotations = annotations[-Model$freq_idx$rem,]
  }
  if (from == "species") {
    Model$annotations = annotations
  } else{
    Model[[from]] = annotations
  }
  return(Model)
}

match_db = function(annotations,
                    db,
                    from = NULL,
                    by_Synonyms = F,
                    Model = NULL) {
  if (is.null(from)) {
    from = colnames(annotations)
  }
  annotations=annotations[,from %in% colnames(db)]
  from = from[from %in% colnames(db)]
  annotations = data.frame(lapply(annotations, toupper))
  db = data.frame(lapply(db, toupper))
  if (is.null(Model)) {
    annotations$idx = 1:length(annotations[, 1])
  } else{
    if (length(Model$freq_idx$rem) > 0) {
      annotations$idx = data.frame(1:nrow(Model$header))[-Model$freq_idx$rem,]
    }
  }
  annotations$ann_idx = 1:length(annotations[, 1])
  annotations$db_idx = NA
  search_subset = annotations
  pb = txtProgressBar(style = 3)
  index = 0
  for (search_item in from) {
    index = index + 1
    setTxtProgressBar(pb, index / length(from))
    search_subset = annotations[!is.na(annotations[search_item]),]
    search_subset = search_subset[is.na(search_subset$db_idx),]
    if (nrow(search_subset) > 0) {
      search_result = lapply(search_subset[[search_item]], function(x) {
        return(which(x == db[[search_item]]))
      })
      for (i in 1:length(search_result)) {
        if (length(search_result[[i]]) > 0) {
          annotations$db_idx[search_subset$ann_idx[i]] = search_result[[i]][1]
        }
      }
    }
  }
  close(pb)
  if (by_Synonyms) {
    search_subset = annotations[!is.na(annotations[search_item]),]
    search_subset = search_subset[is.na(annotations$db_idx),]
    target_list = lapply(db$Synonyms, function(x) {
      return(strsplit(x, ";"))
    })
    pb = txtProgressBar(style = 3)
    search_result = lapply(1:nrow(search_subset), function(i) {
      setTxtProgressBar(pb, i / nrow(search_subset))
      x = search_subset[["Compound.name"]][i]
      for (i in 1:length(target_list)) {
        if (x %in% target_list[[i]][[1]]) {
          return(i)
        }
      }
    })
    close(pb)
    if (length(search_result) > 0) {
      for (i in 1:length(search_result)) {
        if (length(search_result[[i]]) > 0) {
          annotations$db_idx[search_subset$ann_idx[[i]]] = search_result[[i]]
        }
      }
    }
  }
  db_object = db[annotations$db_idx,]
  db_object$idx = annotations$idx
  db_object$db_idx = annotations$db_idx
  rownames(db_object) = 1:length(rownames(db_object))
  return(db_object)
}
database_extension = function(ori_db, ext_db, db_name, by_Synonyms = F) {
  if (!db_name %in% colnames(ori_db)) {
    ori_db[[db_name]] = NA
  }
  for (str in colnames(ori_db)[which(!colnames(ori_db) %in% colnames(ext_db))]) {
    ext_db[str] = NA
  }
  for (str in colnames(ext_db)[which(!colnames(ext_db) %in% colnames(ori_db))]) {
    ext_db = dplyr::select(ext_db,-str)
  }
  ext = match_db(ext_db, ori_db, by_Synonyms = by_Synonyms)
  update_db = ext_db[!is.na(ext$db_idx),]
  update_db$db_idx = ext$db_idx[!is.na(ext$db_idx)]
  ori_db[db_name][update_db$db_idx,] = update_db[db_name]
  add_db = lapply(1:length(ext_db[[1]]), function(x) {
    if (is.na(ext_db[x, db_name])) {
      return(F)
    }
    return(ext_db[x, db_name] %in% ori_db[, db_name])
  })
  add_db = do.call(rbind, add_db)
  add_db = ext_db[!add_db,]
  ret_db = rbind(ori_db, add_db)
  return(ret_db)
}
fill_db_object = function(annotations, db_object) {
  for (i in 1:length(db_object$db_idx)) {
    if (is.na(db_object$db_idx[i])) {
      for (x in colnames(annotations)) {
        if (x %in% colnames(db_object)) {
          db_object[i, which(colnames(db_object) %in% x)] = annotations[i, which(colnames(annotations) %in%
                                                                                   x)]
        }
      }
    }
  }
  return(db_object)
}
Z_score = function(line) {
  return(data.frame((line - mean(line)) / var(line)))
}
get_enzyme = function(Model, db_object,by=c('Uniprot','vmhDB')) {
  by=match.arg(by)
  if(by=="Uniprot"){
    if ("Uniprot.ID" %in% colnames(Model$annotations)) {
      return(db_object[!is.na(Model$annotations$Uniprot.ID),])
    }
    return(db_object[!is.na(db_object$Uniprot.ID),])
  }
}
get_ref_idx = function(meta_data,
                       Model,
                       db,
                       by_Synonyms = F,
                       by_Formula = F) {
  meta_db_obj = match_db(meta_data, db, by_Synonyms = by_Synonyms)
  meta_db_obj = fill_db_object(meta_data, meta_db_obj)
  meta_db_obj = dplyr::select(meta_db_obj,-"idx")
  colnames(meta_db_obj) = gsub("db_idx", "Lab_ID", colnames(meta_db_obj))
  db_object = Model$db_object
  if(!is.null(db_object$idx)){
    db_object = dplyr::select(db_object,-"idx")
  }
  colnames(db_object) = gsub("db_idx", "Lab_ID", colnames(db_object))
  if (!by_Formula) {
    if ("Formula" %in% colnames(db_object)) {
      db_object = dplyr::select(db_object,-"Formula")
    }
    if ("Formula" %in% colnames(meta_db_obj)) {
      meta_db_obj = dplyr::select(meta_db_obj,-"Formula")
    }
  }
  meta_db_obj$level=meta_data$Level
  ref_idx = mutiple_match(db_object,meta_db_obj)
  return(ref_idx)
}
#' Enzyme scoring.
#'
#' formula = "F1" (default, and the definition used throughout the manuscript)
#'   S(E) = sum_i z_i^2 / D_iE^2, i.e. each standardized change is squared
#'   individually and then weighted by the inverse squared network distance.
#' formula = "F0" reproduces the earlier behaviour (sum of the signed weighted
#'   changes, squared afterwards) and is kept only so that previously published
#'   result files can be regenerated.
get_local_score = function(Model,abs=F,do_zscore=T,formula=c("F1","F0")) {
  formula = match.arg(formula)
  identified_meta = Model$ref
  if(abs){
    identified_meta$level=abs(identified_meta$level)
  }
  if(do_zscore){
    identified_meta$level = Z_score(identified_meta$level)
  }
  enzyme = Model$enzyme
  dis_mat = Model$dis_mat
  local_score = lapply(1:nrow(enzyme), function(x) {
    D = dis_mat[, enzyme$idx[x]]
    D = D[as.numeric(identified_meta$idx)]
    D = D[!D == 0]
    w = 1 / (D * D)
    if (formula == "F1") {
      return(sum(w * identified_meta$level^2))
    }
    score = sum(w * identified_meta$level / sum(w))
    return(score^2)
  })
  local_score = data.frame(score = do.call(rbind, local_score))
  local_score[is.na(local_score)] = -10000
  rownames(local_score)=Model$enzyme$Compound.name
  return(local_score)
}
get_p = function(Model,repeat_time,p=NULL,n_cores=10,do_zscore=T,
                 formula=c("F1","F0"),correction=c("add1","old")) {
  formula = match.arg(formula); correction = match.arg(correction)
  pb = txtProgressBar(style = 3)
  if(is.null(p)){
    p = data.frame(p_value = rep(0, nrow(Model$local_score)))
  }
  options(future.globals.maxSize= 891289600*8)
  plan(multisession, workers=n_cores)
  handlers(handler_txtprogressbar(enable = TRUE))
  with_progress({
    p = progressor(along = 1:repeat_time)
    object = future_lapply(
      X = 1:repeat_time,
      FUN = function(i) {
        p()
        setTxtProgressBar(pb, i / repeat_time)
        Model$dis_mat = shuffle_matrix(Model$distance_martix)
        return(
          get_local_score(Model,do_zscore = do_zscore,formula = formula)
        )
      }
    )
  })
  plan(sequential)
  object=do.call(cbind,object)
  temp=lapply(1:ncol(object),function(x){
    return(object[x]>Model$local_score)
  })
  temp=do.call(cbind,temp)
  k=rowSums(temp)
  # correction = "add1" (default): (k+1)/(B+1), the finite-permutation estimator
  # used for every p-value in the manuscript; "old" reproduces max(k,1)/B.
  if (correction == "add1") {
    return((k + 1) / (repeat_time + 1))
  }
  p=lapply(1:length(k),function(x){
    return(max(k[x],1))
  })
  p=do.call(rbind,p)
  return(p / repeat_time)
}
shuffle_matrix = function(matrix) {
  row = nrow(matrix)
  mat2vec = as.matrix(matrix)
  vec2mat = sample(mat2vec, length(mat2vec))
  dim(vec2mat) = c(row, row)
  vec2mat = data.frame(vec2mat)
  return(vec2mat)
}
scatter_plot = function(Model,test_time,
                        p_limit = 0.5,
                        first_n = 5,
                        sign=NULL,
                        save_path=getwd()) {
  plot_data = data.frame(
    local_score = (Model$local_score),
    p_value = Model$p,
    name = Model$enzyme$Compound.name
  )
  plot_data = plot_data[plot_data$local_score >= -1, ]
  plot_data = plot_data[plot_data$p_value <= p_limit, ]
  m1=signif(max(plot_data$local_score), 3)
  m2=signif(min(plot_data$local_score), 3)
  len=(m2-m1)/4
  breaks=c(0,m2+len,m2+2*len,m2+3*len,m1)
  g=ggplot(data = plot_data,)  + guides(fill = guide_legend(byrow = TRUE))+
    geom_point(aes(
      x = local_score,
      y = -log10(p_value),
      color = -log10(p_value),
      size = local_score
    ))+scale_x_continuous(breaks = breaks)+
    scale_size_continuous(range = c(2, 12),guide='none') + 
    scale_color_gradientn(
      values = seq(0, 1, 1 / 7),
      colors = c(
        "#910685",
        "#066CB7",
        "#4CCBC0",
        "#049643",
        "#F7C002",
        "#F39803",
        "#E70015"
      )
    ) +
    theme_classic() +
    xlab("log(Local score)")+
    ylab("-log(P-value)")+
    theme(
      legend.key.height = unit(30, "pt"),
      legend.key.size = unit(35, "pt"),
      legend.key = element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA),
      legend.box.background = element_rect(fill=rgb(1,0,0,alpha=0),colour=NA),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 25),
      legend.margin=unit(0.1,"inches"),
      legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA),
      axis.title = element_text(size = 35),
      axis.text = element_text(size = 30),
      axis.line = element_line(size = 2, colour = "grey30"),
      plot.title = element_text(size=40),
      axis.text.x = element_text(size = 0)
    )+labs(title = paste("Permutation time=",as.character(test_time)))
  g
  if(!is.null(sign)){
    plot_data=plot_data[plot_data$name%in%sign,]
    plot_data[order(plot_data$p_value, decreasing = F)[1:first_n], ]
    g=g+
      geom_text_repel(
        data = plot_data[order(plot_data$p_value, decreasing = F)[1:first_n], ],
        aes(score,-log10(p_value),
            label = name),
        hjust = 3,
        force = 30,
        max.overlaps=100,
        size=10
      )
  }
  else if(first_n>0){
    plot_data=plot_data[order(plot_data$local_score, decreasing = T),]
    plot_data[order(plot_data$p_value, decreasing = F)[1:first_n], ]
    g=g+
      geom_text_repel(
        data = plot_data[order(plot_data$p_value, decreasing = F)[1:first_n], ],
        aes(local_score,-log10(p_value),
            label = name),
        hjust = 2,
        #vjust=-(1:5)*3,
        force = 30,
        max.overlaps=100,
        size=10
      )
  }
  g
  ggsave(paste0(save_path,"\\scatter_plot.pdf"),width=16.9,height = 10.4)
  #dev.off()
}


mutiple_match=function(ref_db,met_db){
  ref_db$idx=rownames(ref_db)
  met_db=met_db[!is.na(met_db$Lab_ID),]
  ref_db=ref_db[!is.na(ref_db$Lab_ID),]
  ret=data.frame()
  pb = txtProgressBar(style = 3)
  for (i in 1:nrow(met_db)) {
    setTxtProgressBar(pb, i / nrow(met_db))
    temp=ref_db[ref_db$Lab_ID==met_db$Lab_ID[i],]
    if(nrow(temp)>0){
      temp$level=met_db$level[i]
      ret=rbind(ret,temp)
    }
  }
  close(pb)
  return(ret)
}

build_model_from_reactions=function(Model){
  Model=get_header(Model)
  model=Model$model@model
  header=data.frame()
  for (i in 1:length(model@reactions)) {
    react=model@reactions[[i]]
    cvTerms=react@cvTerms
    cvs=c()
    for (cv in cvTerms) {
      for (j in cv@resources) {
        aux=strsplit(j,'/')[[1]]
        if(aux[length(aux)-1]=='ec-code'){
          cvs=c(cvs,aux[length(aux)])
          model@reactions[[i]]@modifiers[[length(model@reactions[[i]]@modifiers)+1]]=1
          names(model@reactions[[i]]@modifiers)[length(model@reactions[[i]]@modifiers)]=aux[length(aux)]
        }
      }
    }
    if(!is.null(cvs)){
      header=rbind(header,data.frame(com_name=cvs,compartment="ud"))
    }
  }
  header=data.frame(com_name=unique(header$com_name),compartment='ud')
  rownames(header)=header$com_name
  Model$enzyme=data.frame(idx=(nrow(Model$header)+1):(length(header$com_name)+nrow(Model$header)),
                          Compound.name=header$com_name)
  Model$header=rbind(Model$header,header)
  Model$model@model=model
  return(Model)
}

metabolism_annotation=function(path,polar,database,filter=NA,use_MS2=F,min_fraction = 0,column = "rp",snthresh = 5,noise=500){
  setwd(path)
  massprocesser::process_data(
    path = "./",
    polarity = polar,
    ppm = 25,
    peakwidth = c(5, 30),
    snthresh = snthresh,
    noise = noise,
    threads = 5,
    output_tic = F,
    output_bpc = F,
    output_rt_correction_plot = F,
    min_fraction = min_fraction,
    fill_peaks = FALSE
  )
  setwd(paste0(path,'\\Result'))
  raw_data = read.csv("peak_table_for_cleaning.csv",row.names = 1,header = T)
  expression_data = data.frame(raw_data[,-1:-2])
  variable_info = raw_data[,1:2]%>%
    mutate(variable_id=rownames(raw_data))
  variable_info = variable_info[,c(3,1,2)]
  row.names(variable_info) =1:length(row.names(raw_data))
  sample_info = data.frame(sample_id = colnames(expression_data),
                            class= "Subject",
                            group=gsub("[0-9]", "", colnames(expression_data)))
  object =
    create_mass_dataset(
      expression_data = expression_data,
      sample_info = sample_info,
      variable_info = variable_info
    )
  if(use_MS2){
    object =
    mutate_ms2(
      object = object,
      column = column,
      polarity = polar,
      ms1.ms2.match.mz.tol = 15,
      ms1.ms2.match.rt.tol = 30,
      path = path
    )
  }
  object %>%
    activate_mass_dataset(what = "sample_info") %>%
    dplyr::count(group)
  object = 
    impute_mv(object = object, method = "knn")
  used = apply(object@expression_data, 1, function(row){
    !any(row<500)
  })
  object@expression_data = object@expression_data[used,]
  object2 =
    normalize_data(object, method = "median")
  write.csv(object2@expression_data,"knn_expression_data.csv")
  data_base = database
  data_base@spectra.info$RT = NA
  object = object2
  object1 =
    annotate_metabolites_mass_dataset(object = object,
                                      database = data_base,
                                      ms1.match.ppm = 15,
                                      threads =3,
                                      polarity=polar,
                                      candidate.num = 5)
  annotations = extract_annotation_table(object1) 
  write.xlsx(annotations,"All_annotations.xlsx",row.names = F,fileEncoding='UTF-8')
  filted_annotations=data.frame()
  if(!is.na(filter)){
    for (i in filter) {
      filted_annotations=rbind(filted_annotations,annotations[annotations$Adduct==i,])
    }
  }else{
    app1=annotations
  }
  filted_annotations=filted_annotations[,!colnames(filted_annotations)=="Database"]
  write.xlsx(filted_annotations,"filted.xlsx",row.names = F,fileEncoding='UTF-8')
}

ec2gen=function(x){
  a=lapply(x,function(x){
    ret=c()
    tryCatch({
      ret=KEGGREST::keggGet(x)
    }, error = function(e){
      print(x)
    })
    return(ret)
  })
  gene_kidney=lapply(a,function(x){
    G=x[[1]]$GENES[1]
    if(!is.null(G)){
      if(substr(G,1,3)=="HSA"){
        left=c()
        right=c()
        for (i in 5:nchar(G)) {
          if(substr(G,i,i)=="("){
            left=c(left,i)
          }else if(substr(G,i,i)==")"){
            right=c(right,i)
          }
        }
        ret=c()
        for (i in 1:length(left)) {
          ret=c(ret,substr(G,left[i]+1,right[i]-1))
        }
        return(ret)
      }
    }
    return(NA)
  }
  )
}
getresult=function(Model){
  return(data.frame(score=Model$local_score,p=Model$p,enzyme=Model$enzyme$Compound.name))
}
set_label = function(query_id, database) {
  return_list = list()
  for (i in seq_along(database)) {
    x = database[[i]] 
    temp = match(query_id, x)
    temp[!is.na(temp)] = 1
    temp[is.na(temp)] = 0
    return_list[[i]] = temp
  }
  return_result = do.call(cbind, return_list) %>% as.data.frame()
  rownames(return_result) = query_id
  return(return_result)
}
enrichment=function(Model,p_lim=0.05,method = c("hypergeometric", "fisher"), 
                    p_adjust_method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                    plot_dir=getwd(),plot_axis="p_value",pdf_name="Enrich_informations"){
  match.arg(db)
  match.arg(method)
  match.arg(p_adjust_method)
  setwd(plot_dir)
  pathway_database=kegg_hsa_pathway
  result=getresult(Model)
  result=result[result$p<p_lim,]
  result=result[order(result$p),]
  query_id=result$enzyme
  if (pathway_database@database_info$source != "KEGG") {
    stop("pathway_database must from KEGG.\n")
  }
  remain_idx = pathway_database@compound_list %>%
    purrr::map(function(x) {
      nrow(x)
    }) %>%
    unlist() %>%
    (function(x) x > 0)%>% 
    which()
  pathway_database = filter_pathway(object = pathway_database, 
                                    remain_idx = remain_idx)
  for (i in 1:length(pathway_database@pathway_id)) {
    pathway_database@compound_list[[i]] = pathway_database@compound_list[[i]] %>% 
      dplyr::filter(!is.na(KEGG.ID)) %>% dplyr::filter(KEGG.ID != "")
  }
  remain_idx = pathway_database@compound_list %>%
    purrr::map(function(x) {
      nrow(x)
    }) %>%
    unlist() %>%
    (function(x) x > 0) %>% 
    which()
  pathway_database = filter_pathway(object = pathway_database, 
                                    remain_idx = remain_idx)
  if (length(pathway_database@pathway_id) == 0) {
    stop("Pathway database length is zero.\n")
  }else {
    message(crayon::green(length(pathway_database@describtion), 
                          "pathways."))
  }
  database = pathway_database@gene_list
  names(database) = pathway_database@pathway_id
  database = database %>% purrr::map(function(x) {
    temp=unique(do.call(c,strsplit(str_extract(x$Gene.name, "(?<=EC:).*?(?=]|$)")," ")))
    temp=temp[!is.na(temp)]
  })
  remain_idx=c()
  for (i in 1:length(database)) {
    if(length(database[[i]])>0){
      remain_idx=c(remain_idx,i)
    }
  }
  pathway_database =
    filter_pathway(object = pathway_database,
                   remain_idx = remain_idx)
  database=database[remain_idx]
  for (idx in 1:length(pathway_database@describtion)) {
    if (is.null(pathway_database@describtion[[idx]])) {
      pathway_database@describtion[[idx]] = NA
    }
    pathway_database@describtion[[idx]] = paste(pathway_database@describtion[[idx]], 
                                                collapse = "{}")
  }
  pathway_info = data.frame(pathway_id = pathway_database@pathway_id, 
                             pathway_name = pathway_database@pathway_name, describtion = unlist(pathway_database@describtion), 
                             pathway_class = unlist(pathway_database@pathway_class), 
                             stringsAsFactors = FALSE)
  all_id = database %>% unlist() %>% unique() %>% unname() %>% 
    as.character()
  query_id = query_id[which(is.element(query_id, all_id))]
  query_id = unique(query_id[!is.na(query_id)])
  if (length(query_id) == 0) {
    return(NULL)
  }
  sig_id = as.character(query_id)
  num_all = length(all_id)
  num_sig = length(sig_id)
  all_matrix = set_label(query_id = all_id, database = database)
  all_matrix2 = all_matrix[, colSums(all_matrix) != 0, drop = FALSE]
  if (ncol(all_matrix2) < 2) {
    return(NULL)
  }
  sig_matrix = set_label(query_id = sig_id, database = database)
  sig_matrix = sig_matrix[, colSums(all_matrix) != 0, drop = FALSE]
  p_value = furrr::future_map(.x = 1:ncol(all_matrix2), .f = function(i) {
    a1 = sum(sig_matrix[, i])
    a2 = sum(all_matrix2[, i]) - sum(sig_matrix[, i])
    a3 = length(sig_id) - a1
    a4 = (length(all_id) - length(sig_id)) - a2
    tab = t(matrix(c(a1, a2, a3, a4), 2))
    if (method == "hypergeometric") {
      phyper(q = a1 - 1, m = sum(all_matrix2[, i]), n = num_all - 
               sum(all_matrix2[, i]), k = num_sig, lower.tail = FALSE)
    }else {
      check = tryCatch({
        resfish = fisher.test(tab, alternative = "greater")
      }, error = function(e) {
        NA
      })
      if (!is(object = check, class2 = "htest")) {
        return(1)
      }
      else {
        resfish = fisher.test(tab, alternative = "greater")
        return(resfish$p.value)
      }
    }
  }, .progress = TRUE) %>% unlist()
  p_value_adjust = p.adjust(p_value, method = p_adjust_method)
  result = data.frame(pathway_info, p_value, p_value_adjust)
  result$all_id = furrr::future_map(.x = result$pathway_id,.f = function(x) {
    paste(database[[x]], collapse = ";")
  }) %>% unlist()
  result$all_number = furrr::future_map(.x = result$pathway_id,.f = function(x) {
    length(database[[x]])
  }) %>% unlist()
  result$mapped_id = furrr::future_map(.x = result$pathway_id,.f = function(x) {
    paste(query_id[query_id %in% database[[x]]], collapse = ";")
  }) %>% unlist()
  result$mapped_number = furrr::future_map(.x = result$pathway_id,.f = function(x) {
    sum(query_id %in% database[[x]])
  }) %>% unlist()
  result$mapped_percentage = furrr::future_map(.x = result$pathway_id,.f = function(x) {
    sum(query_id %in% database[[x]]) * 100/length(database[[x]])
  }) %>% unlist()
  result = result %>% dplyr::arrange(p_value_adjust)
  return_result = new(Class = "enrich_result", pathway_database = pathway_database@database_info$source, 
                      pathway_version = '1.0.8', result = result, 
                      parameter = new(Class = "tidymass_parameter", pacakge_name = "metpath", 
                                      function_name = "enrich_kegg()",
                                      parameter = list(query_id = query_id, 
                                                       query_type = "gene", id_type = "KEGG", 
                                                       pathway_database = paste(pathway_database@database_info$source, 
                                                                                pathway_database@database_info$version, sep = ","), 
                                                       p_cutoff = p_lim, p_adjust_method = p_adjust_method, 
                                                       method = method, threads = 5), time = Sys.time()))
  plot_data = return_result@result %>%
    mutate(
      enrichment_ratio = mapped_number / all_number
    )
  if (plot_axis == "p_value") {
    plot_data = plot_data %>%
      filter(p_value < 0.05) %>%
      mutate(neg_log_p = -log10(p_value))
  } else if (plot_axis == "p_value_adjust") {
    if (!"p_value_adjust" %in% names(plot_data)) {
      stop("Error: When plot_axis is 'p_value_adjust', the 'p_value_adjust' column was not found in the 'result' dataframe.")
    }
    plot_data = plot_data %>%
      filter(p_value_adjust < 0.05) %>%
      mutate(neg_log_p = -log10(p_value_adjust))
  } else {
    stop("Invalid 'plot_axis' value. It must be 'p_value' or 'p_value_adjust'.")
  }
  
  if (nrow(plot_data) == 0) {
    message(paste0("No pathway enriched"))
  }
  
  plot_data = plot_data %>%
    mutate(pathway_name = fct_reorder(pathway_name, enrichment_ratio))
  lollipop_plot = ggplot(plot_data, aes(x = enrichment_ratio, y = pathway_name)) +
    geom_segment(aes(x = 0, xend = enrichment_ratio, y = pathway_name, yend = pathway_name),
                 color = "grey", linewidth =1) +
    geom_point(aes(size = neg_log_p), shape = 21, fill = "steelblue", stroke = 0) +
    scale_size_continuous(range = c(3, 10), name = paste0("-log10(", gsub("_", " ", plot_axis), ")")) +
    labs(
      x = "Enrichment ratio",
      y = "pathways",
      title = "Enrichment results"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(hjust = 1,size=12),
      axis.text.y = element_text(size=12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    );lollipop_plot
  pdf(paste0("Enrich_informations",".pdf"), width = 8, height = 10,family = "sans")
  enrich_scatter_plot(object=return_result,y_axis=plot_axis)
  lollipop_plot
  dev.off()
  return(return_result)
}
