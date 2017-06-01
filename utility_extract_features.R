library(biomaRt)

gene_ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
			host="grch37.ensembl.org", 
			path="/biomart/martservice" ,
			dataset="hsapiens_gene_ensembl")

funcgene_motif_ensembl = useMart(
			biomart="ENSEMBL_MART_FUNCGEN",
			host="grch37.ensembl.org",
			path="/biomart/martservice",
			dataset="hsapiens_motif_feature")

funcgene_regulatory_ensembl = useMart(biomart="ENSEMBL_MART_FUNCGEN",
				host="grch37.ensembl.org",
				path="/biomart/martservice",
				dataset="hsapiens_regulatory_feature")

hum_ext_feature = useMart(biomart="ENSEMBL_MART_FUNCGEN",
			host="grch37.ensembl.org",
			path="/biomart/martservice",
			dataset="hsapiens_external_feature")

strvariant_ensembl = useMart(biomart="ENSEMBL_MART_SNP",
				host="grch37.ensembl.org",
				path="/biomart/martservice",
				dataset="hsapiens_structvar")

#snp_ensembl = useMart(biomart="ENSEMBL_MART_SNP", 
#			host="grch37.ensembl.org",
#			path="/biomart/martservice",
#			dataset="hsapiens_snp")

get_gene_info <- function(label,chr,downstream_flank=0,upstream_flank=0){
  tbl_info <- NA
  switch(label,
	 exon={
	       tbl_info <- getBM(attributes=c('chromosome_name',
               			        'external_gene_name',
                       			'ensembl_transcript_id',
                       			'exon_chrom_start',
                       			'exon_chrom_end','strand'),
               			filters=c('chromosome_name'),
               			value=list(chr),mart=gene_ensembl)
  		
		names(tbl_info)[4:5] <- c("coord_start","coord_end")

	 },
	 intron={
		tmp_tbl <- getBM(attributes=c('chromosome_name',
                                        'external_gene_name',
                                        'ensembl_transcript_id',
                                        'exon_chrom_start',
                                        'exon_chrom_end','strand'),
                                filters=c('chromosome_name'),
                                value=list(chr),mart=gene_ensembl)
		
		for(transcrpt in unique(tmp_tbl$ensembl_transcript_id)){
		    
		   strand=unique(tmp_tbl[which(tmp_tbl$ensembl_transcript_id == transcrpt),"strand"])
		   t_tbl <- tmp_tbl[which(tmp_tbl$ensembl_transcript_id == transcrpt),]
		   t_tbl <- t_tbl[order(t_tbl$exon_chrom_start,t_tbl$exon_chrom_start),]
		   if(dim(t_tbl)[1]  < 2) next;
		   
		   exon_tbl <- cbind(t_tbl$exon_chrom_start,t_tbl$exon_chrom_end)
		   exon_tbl <- exon_tbl[order(exon_tbl[,1]),]
		   exon_vec <- as.vector(unlist(exon_tbl))
		   exon_vec <- exon_vec[-c(1,length(exon_vec))]
		   intron_tbl <- data.frame(matrix(exon_vec,ncol=2))	
		   intron_tbl <- data.frame(
					   chromosome_name=chr,
					   ensembl_transcript_id=transcrpt,	  	
					   coord_start =intron_tbl[,2]+1,
					   coord_end=intron_tbl[,1]-1,
					   strand=strand
					)
				
		  tbl_info <- rbind(tbl_info,intron_tbl) 
		}
		rm(tmp_tbl)	

	 },
         five_prime_utr={
		 tbl_info <- getBM(attributes=c('chromosome_name',
                                	'ensembl_gene_id',
                                	'ensembl_transcript_id',
                               		'gene_biotype',
                                	'5_utr_start',
                                	'5_utr_end','strand'),
               			filters=c('chromosome_name'),
                		value=list(chr),mart=gene_ensembl)
  		  names(tbl_info)[5:6] <- c("coord_start","coord_end")
	 },
	 three_prime_utr={
                 tbl_info <- getBM(attributes=c('chromosome_name',
                                        'ensembl_gene_id',
                                        'ensembl_transcript_id',
                                        'gene_biotype',
                                        '3_utr_start',
                                        '3_utr_end','strand'),
                                filters=c('chromosome_name'),
                                value=list(chr),mart=gene_ensembl)
                  names(tbl_info)[5:6] <- c("coord_start","coord_end")
         },
	 mirna={
		tbl_info <- getBM(attributes=c('chromosome_name',
                               'mirbase_id',
                               #'ensembl_gene_id',
                                'start_position',
                                'end_position',
                                'strand'),
               		filters=c('chromosome_name'),
                	value=list(chr),mart=gene_ensembl)
	
		tbl_info <- tbl_info[!(is.na(tbl_info$mirbase_id)|tbl_info$mirbase_id == ""),]
		names(tbl_info)[4:5] <- c("coord_start","coord_end")
	 },
	 TF_binding_motif={
                tbl_info <- getBM(attributes=c( 'binding_matrix_id',
                                                'chromosome_name',
                                                'chromosome_start',
                                                'chromosome_end',
						'chromosome_strand',
                                                'so_accession',
                                                'so_name'
                                                ),
                                filters=c('chromosome_name'),
                                value=list(chr),mart=funcgene_motif_ensembl)
                names(tbl_info)[3:5] <- c("coord_start","coord_end","strand")
        },
	StructuralVariant={
	  	tbl_info <- getBM(attributes=c('dgva_study_accession',
                                   		'sv_accession',
                                   		'sv_variant_type',
						'allele_string',
						'structural_variation_clinical_significance',
                                   		'chr_name',
                                   		'chrom_start',
                                   		'chrom_start',
				   		'seq_region_strand'),
				filters=c('chr_name'),
                                value=list(chr),mart=strvariant_ensembl)

		names(tbl_info)[6:9] <- c("chromosome_name","coord_start","coord_end","strand")	
	},      
	{
		stop(cat("Function Not defined for label:",label,"\n",
                "Labels defined for get_info:","\n",
                "exon","\n",
                "intron","\n",
                "five_prime_utr","\n",
                "three_prime_utr","\n",
                "mirna","\n",
		"TF_binding_motif","\n",
		"StructuralVariant","\n"))
	}
      )
   tbl_info=tbl_info[order(tbl_info$coord_start ),]
   select_pos=tbl_info$strand==1
   tbl_info$coord_start[select_pos] = tbl_info$coord_start[select_pos] - downstream_flank
   tbl_info$coord_end[select_pos]  = tbl_info$coord_end[select_pos] + upstream_flank

   tbl_info$coord_start[!select_pos] = tbl_info$coord_start[!select_pos] + downstream_flank
   tbl_info$coord_end[!select_pos]  = tbl_info$coord_end[!select_pos] - upstream_flank

   tbl_info=tbl_info[order(tbl_info$coord_start ),]
   return(tbl_info)
}

get_regulatory_info <- function(label,chr){
  tbl_info <- NA
  switch(label,
	 enhancer={
		tbl_info <- getBM(attributes=c('chromosome_name',
                                 		'chromosome_start',
                                 		'chromosome_end',
						'so_accession',
                                 		'so_name'),
				filters=c('chromosome_name'),
				value=list(chr),mart=hum_ext_feature)
		tbl_info <- tbl_info[which(tbl_info$so_name == "enhancer"),]
		names(tbl_info)[2:3] <- c("coord_start","coord_end")

	},
	#Transcription Start Site
	transcription_start_site={
		tbl_info <- getBM(attributes=c('chromosome_name',
                                                'chromosome_start',
                                                'chromosome_end',
                                                'so_accession',
                                                'so_name'),
                                filters=c('chromosome_name'),
                                value=list(chr),mart=hum_ext_feature)
                tbl_info <- tbl_info[which(tbl_info$so_name == "transcription_start_site"),]
                names(tbl_info)[2:3] <- c("coord_start","coord_end")

	},
	promoter={
		tbl_info <- getBM(attributes=c('chromosome_name',
						'bound_seq_region_start',
						'bound_seq_region_end',		
                                                'chromosome_start',
                                                'chromosome_end',
                                                'so_accession',
                                                'so_name'),
                                filters=c('chromosome_name'),
                                value=list(chr),mart=funcgene_regulatory_ensembl)
		tbl_info <- tbl_info[which(tbl_info$so_name == "promoter"),]
		names(tbl_info)[4:5] <- c("coord_start","coord_end")
	},
	TF_binding_site={
		tbl_info <- getBM(attributes=c('chromosome_name',
						'bound_seq_region_start',
                                                'bound_seq_region_end',
                                                'chromosome_start',
                                                'chromosome_end',
                                                'so_accession',
                                                'so_name'),
                                filters=c('chromosome_name'),
                                value=list(chr),mart=funcgene_regulatory_ensembl)
                tbl_info <- tbl_info[which(tbl_info$so_name == "TF_binding_site"),]
                names(tbl_info)[4:5] <- c("coord_start","coord_end")
	},
	open_chromatin_region={
		 tbl_info <- getBM(attributes=c('chromosome_name',
						'bound_seq_region_start',
                                                'bound_seq_region_end',	
                                                'chromosome_start',
                                                'chromosome_end',
                                                'so_accession',
                                                'so_name'),
                                filters=c('chromosome_name'),
                                value=list(chr),mart=funcgene_regulatory_ensembl)
                tbl_info <- tbl_info[which(tbl_info$so_name == "open_chromatin_region"),]
                names(tbl_info)[4:5] <- c("coord_start","coord_end")

	},
        {
	    stop(cat("Function Not defined for label:",label,"\n",
		"Labels defined for get_info:","\n",
		"enhancer","\n",
		"transcription_start_site","\n",
		"promoter","\n",
		"TF_binding_site","\n",
		"open_chromatin_region","\n"		
			))
	 }     
        ) 
   return(tbl_info)
}

