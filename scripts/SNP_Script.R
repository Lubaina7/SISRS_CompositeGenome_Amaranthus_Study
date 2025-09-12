library(Rboretum)
sourceRboretum()

# species_signal <- getAlignmentSignal('data/raw_data/Species pooled/alignment_m0_nogap.phylip-relaxed') %>% mutate(Loc = read_lines('data/raw_data/Species pooled/alignment_locs_m0_nogap.txt'))

# species_signal %>%
  # write_tsv('data/output/Species_Signal.tsv',quote = "none")

# ind_signal <- getAlignmentSignal('data/raw_data/Sample specific/alignment_pi_singletons_m17_nogap.phylip-relaxed') %>% mutate(Loc = read_lines('data/raw_data/Sample specific/alignment_pi_singletons_locs_m17_nogap.txt'))

# ind_signal %>%
# write_tsv('data/output/Individual_Signal.tsv',quote = "none")

species_signal <- read_tsv('data/output/Species_Signal.tsv')

ind_signal <- read_tsv('data/output/Individual_Signal.tsv') %>% 
  mutate(Loc = read_lines('data/raw_data/Sample specific/alignment_pi_singletons_locs_m17_nogap.txt'))

species_id_df <- read_csv('data/ID_Conversion.csv')

contig_lengths <- read_tsv('data/raw_data/Sample specific/contigs_SeqLength.tsv') %>% `colnames<-`(c("Contig","Length"))

# Species lists
atub_samples <- species_id_df %>% filter(Species_Pool == "Atub") %>% pull(Sample_ID) %>% sort()
awat_samples <- species_id_df %>% filter(Species_Pool == "Awat") %>% pull(Sample_ID) %>% sort()
aret_samples <- species_id_df %>% filter(Species_Pool == "Aret") %>% pull(Sample_ID) %>% sort()
apow_samples <- species_id_df %>% filter(Species_Pool == "Apow") %>% pull(Sample_ID) %>% sort()
ahyp_samples <- species_id_df %>% filter(Species_Pool == "Ahyp") %>% pull(Sample_ID) %>% sort()
aare_samples <- species_id_df %>% filter(Species_Pool == "Aare") %>% pull(Sample_ID) %>% sort()
apal_samples <- species_id_df %>% filter(Species_Pool == "Apal") %>% pull(Sample_ID) %>% sort()
atri_samples <- species_id_df %>% filter(Species_Pool == "Atri") %>% pull(Sample_ID) %>% sort()
acau_samples <- species_id_df %>% filter(Species_Pool == "Acau") %>% pull(Sample_ID) %>% sort()

# Species strings
atub_string <- vectorSemi(atub_samples)
awat_string <- vectorSemi(awat_samples)
aret_string <- vectorSemi(aret_samples)
apow_string <- vectorSemi(apow_samples)
ahyp_string <- vectorSemi(ahyp_samples)
aare_string <- vectorSemi(aare_samples)
apal_string <- vectorSemi(apal_samples)
atri_string <- vectorSemi(atri_samples)
acau_string <- vectorSemi(acau_samples)

# Clades of interest


# 28,275
upper_clade_string <- species_id_df %>% filter(Upper_Clade == 1) %>% pull(Sample_ID) %>% sort() %>% vectorSemi()
upper_clade_locs <- ind_signal %>%
  filter(Split_1 == upper_clade_string | Split_2 == upper_clade_string | Split_3 == upper_clade_string | Split_4 == upper_clade_string) %>%
  mutate(SNP_Species = "Upper_Clade")

# 17,990 
tubare_string <- species_id_df %>% filter(TubAre == 1) %>% pull(Sample_ID) %>% sort() %>% vectorSemi()
tubare_clade_locs <- ind_signal %>%
  filter(Split_1 == tubare_string | Split_2 == tubare_string | Split_3 == tubare_string | Split_4 == tubare_string) %>%
  mutate(SNP_Species = "TubAre")

# 47,523
lower_clade_a_string <- species_id_df %>% filter(Lower_Clade_A == 1) %>% pull(Sample_ID) %>% sort() %>% vectorSemi()
lower_clade_a_clade_locs <- ind_signal %>%
  filter(Split_1 == lower_clade_a_string | Split_2 == lower_clade_a_string | Split_3 == lower_clade_a_string | Split_4 == lower_clade_a_string) %>%
  mutate(SNP_Species = "Lower_Clade_A")

# 70,105
lower_clade_b_string <- species_id_df %>% filter(Lower_Clade_B == 1) %>% pull(Sample_ID) %>% sort() %>% vectorSemi()
lower_clade_b_clade_locs <- ind_signal %>%
  filter(Split_1 == lower_clade_b_string | Split_2 == lower_clade_b_string | Split_3 == lower_clade_b_string | Split_4 == lower_clade_b_string) %>%
  mutate(SNP_Species = "Lower_Clade_B")

# 22,031
caud_hypo_string <- species_id_df %>% filter(CaudHypo == 1) %>% pull(Sample_ID) %>% sort() %>% vectorSemi()
caud_hypo_clade_locs <- ind_signal %>%
  filter(Split_1 == caud_hypo_string | Split_2 == caud_hypo_string | Split_3 == caud_hypo_string | Split_4 == caud_hypo_string) %>%
  mutate(SNP_Species = "CaudHypo")

# 9,300
pow_ret_string <- species_id_df %>% filter(PowRet == 1) %>% pull(Sample_ID) %>% sort() %>% vectorSemi()
pow_ret_string2 <- species_id_df %>% filter(PowRet == 1) %>% filter(Sample_ID != "trimmed_22002D-01-03_S7_L001_Nuclear")%>% pull(Sample_ID) %>% sort() %>% vectorSemi()
pow_ret_clade_locs <- ind_signal %>%
  filter(Split_1 == pow_ret_string | Split_2 == pow_ret_string | Split_3 == pow_ret_string | Split_4 == pow_ret_string |
           Split_1 == pow_ret_string2 | Split_2 == pow_ret_string2 | Split_3 == pow_ret_string2 | Split_4 == pow_ret_string2) %>%
  mutate(SNP_Species = "PowRet")

# all_big_locs <- bind_rows(upper_clade_locs,tubare_clade_locs,lower_clade_a_clade_locs,lower_clade_b_clade_locs,caud_hypo_clade_locs,pow_ret_clade_locs) %>%
#   separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE)
# 
# one_major_snp <- all_big_locs %>% 
#   rowwise() %>%
#   filter(
#     !all(atub_samples %in% semiVector(Non_Base_Taxa)) &
#       !all(awat_samples %in% semiVector(Non_Base_Taxa)) &
#       !all(aret_samples %in% semiVector(Non_Base_Taxa)) &
#       !all(apow_samples %in% semiVector(Non_Base_Taxa)) &
#       !all(ahyp_samples %in% semiVector(Non_Base_Taxa)) &
#       !all(aare_samples %in% semiVector(Non_Base_Taxa)) &
#       !all(apal_samples %in% semiVector(Non_Base_Taxa)) &
#       !all(atri_samples %in% semiVector(Non_Base_Taxa)) &
#       !all(acau_samples %in% semiVector(Non_Base_Taxa))) %>%
#   filter(Non_Base_Count < 5) %>%
#   group_by(Contig) %>%
#   count(SNP_Species) %>%
#   pivot_wider(names_from = "SNP_Species",values_from = n,values_fill = 0) %>% 
#   rowwise() %>% 
#   mutate(Major_Min = min(c(Lower_Clade_B,Lower_Clade_A,TubAre,Upper_Clade)),
#          Minor_Min = min(c(PowRet,CaudHypo))) %>% 
#   filter(Major_Min >= 2 & Minor_Min >= 1)
# 


# Species SNPs

# atub_fixed_locs <- ind_signal %>%
#   rowwise() %>%
#   filter(all(atub_samples %in% semiVector(Split_1)) | 
#            all(atub_samples %in% semiVector(Split_2)) |
#            all(atub_samples %in% semiVector(Split_3)) |
#            all(atub_samples %in% semiVector(Split_4))) %>%
#   ungroup()
# 
# aare_fixed_locs <- ind_signal %>%
#   rowwise() %>%
#   filter(all(aare_samples %in% semiVector(Split_1)) | 
#            all(aare_samples %in% semiVector(Split_2)) |
#            all(aare_samples %in% semiVector(Split_3)) |
#            all(aare_samples %in% semiVector(Split_4))) %>%
#   ungroup()
# 
# 
# 
# 
# awat_fixed_locs <- ind_signal %>%
#   rowwise() %>%
#   filter(all(awat_samples %in% semiVector(Split_1)) | 
#            all(awat_samples %in% semiVector(Split_2)) |
#            all(awat_samples %in% semiVector(Split_3)) |
#            all(awat_samples %in% semiVector(Split_4))) %>%
#   ungroup()
# 
# aret_fixed_locs <- ind_signal %>%
#   rowwise() %>%
#   filter(all(aret_samples %in% semiVector(Split_1)) | 
#            all(aret_samples %in% semiVector(Split_2)) |
#            all(aret_samples %in% semiVector(Split_3)) |
#            all(aret_samples %in% semiVector(Split_4))) %>%
#   ungroup()
# 
# apow_fixed_locs <- ind_signal %>%
#   rowwise() %>%
#   filter(all(apow_samples %in% semiVector(Split_1)) | 
#            all(apow_samples %in% semiVector(Split_2)) |
#            all(apow_samples %in% semiVector(Split_3)) |
#            all(apow_samples %in% semiVector(Split_4))) %>%
#   ungroup()
# 
# ahyp_fixed_locs <- ind_signal %>%
#   rowwise() %>%
#   filter(all(ahyp_samples %in% semiVector(Split_1)) | 
#            all(ahyp_samples %in% semiVector(Split_2)) |
#            all(ahyp_samples %in% semiVector(Split_3)) |
#            all(ahyp_samples %in% semiVector(Split_4))) %>%
#   ungroup()
# 
# 
# 
# apal_fixed_locs <- ind_signal %>%
#   rowwise() %>%
#   filter(all(apal_samples %in% semiVector(Split_1)) | 
#            all(apal_samples %in% semiVector(Split_2)) |
#            all(apal_samples %in% semiVector(Split_3)) |
#            all(apal_samples %in% semiVector(Split_4))) %>%
#   ungroup()
# 
# atri_fixed_locs <- ind_signal %>%
#   rowwise() %>%
#   filter(all(atri_samples %in% semiVector(Split_1)) | 
#            all(atri_samples %in% semiVector(Split_2)) |
#            all(atri_samples %in% semiVector(Split_3)) |
#            all(atri_samples %in% semiVector(Split_4))) %>%
#   ungroup()
# 
# acau_fixed_locs <- ind_signal %>%
#   rowwise() %>%
#   filter(all(acau_samples %in% semiVector(Split_1)) | 
#            all(acau_samples %in% semiVector(Split_2)) |
#            all(acau_samples %in% semiVector(Split_3)) |
#            all(acau_samples %in% semiVector(Split_4))) %>%
#   ungroup()
# 
# # Save fixed locs
# atub_fixed_locs %>%
#   rowwise() %>%
#   mutate(Split = ifelse(any(atub_samples %in% semiVector(Split_1)),"Split_1",
#                         ifelse(any(atub_samples %in% semiVector(Split_2)),"Split_2",
#                                ifelse(any(atub_samples %in% semiVector(Split_3)),"Split_3","Split_4")))) %>%
#   ungroup() %>%
#   mutate(Species = "Atub") %>%
#   select(Species,Loc,Split) %>%
#   write_tsv('data/output/Fixed_Locs/Atub.tsv')
# 
# aare_fixed_locs %>%
#   rowwise() %>%
#   mutate(Split = ifelse(any(aare_samples %in% semiVector(Split_1)),"Split_1",
#                         ifelse(any(aare_samples %in% semiVector(Split_2)),"Split_2",
#                                ifelse(any(aare_samples %in% semiVector(Split_3)),"Split_3","Split_4")))) %>%
#   ungroup() %>%
#   mutate(Species = "aare") %>%
#   select(Species,Loc,Split) %>%
#   write_tsv('data/output/Fixed_Locs/Aare.tsv')
# 
# awat_fixed_locs %>%
#   rowwise() %>%
#   mutate(Split = ifelse(any(awat_samples %in% semiVector(Split_1)),"Split_1",
#                         ifelse(any(awat_samples %in% semiVector(Split_2)),"Split_2",
#                                ifelse(any(awat_samples %in% semiVector(Split_3)),"Split_3","Split_4")))) %>%
#   ungroup() %>%
#   mutate(Species = "awat") %>%
#   select(Species,Loc,Split) %>%
#   write_tsv('data/output/Fixed_Locs/Awat.tsv')
# 
# aret_fixed_locs %>%
#   rowwise() %>%
#   mutate(Split = ifelse(any(aret_samples %in% semiVector(Split_1)),"Split_1",
#                         ifelse(any(aret_samples %in% semiVector(Split_2)),"Split_2",
#                                ifelse(any(aret_samples %in% semiVector(Split_3)),"Split_3","Split_4")))) %>%
#   ungroup() %>%
#   mutate(Species = "aret") %>%
#   select(Species,Loc,Split) %>%
#   write_tsv('data/output/Fixed_Locs/Aret.tsv')
# 
# apow_fixed_locs %>%
#   rowwise() %>%
#   mutate(Split = ifelse(any(apow_samples %in% semiVector(Split_1)),"Split_1",
#                         ifelse(any(apow_samples %in% semiVector(Split_2)),"Split_2",
#                                ifelse(any(apow_samples %in% semiVector(Split_3)),"Split_3","Split_4")))) %>%
#   ungroup() %>%
#   mutate(Species = "apow") %>%
#   select(Species,Loc,Split) %>%
#   write_tsv('data/output/Fixed_Locs/Apow.tsv')
# 
# ahyp_fixed_locs %>%
#   rowwise() %>%
#   mutate(Split = ifelse(any(ahyp_samples %in% semiVector(Split_1)),"Split_1",
#                         ifelse(any(ahyp_samples %in% semiVector(Split_2)),"Split_2",
#                                ifelse(any(ahyp_samples %in% semiVector(Split_3)),"Split_3","Split_4")))) %>%
#   ungroup() %>%
#   mutate(Species = "ahyp") %>%
#   select(Species,Loc,Split) %>%
#   write_tsv('data/output/Fixed_Locs/Ahyp.tsv')
# 
# apal_fixed_locs %>%
#   rowwise() %>%
#   mutate(Split = ifelse(any(apal_samples %in% semiVector(Split_1)),"Split_1",
#                         ifelse(any(apal_samples %in% semiVector(Split_2)),"Split_2",
#                                ifelse(any(apal_samples %in% semiVector(Split_3)),"Split_3","Split_4")))) %>%
#   ungroup() %>%
#   mutate(Species = "apal") %>%
#   select(Species,Loc,Split) %>%
#   write_tsv('data/output/Fixed_Locs/Apal.tsv')
# 
# atri_fixed_locs %>%
#   rowwise() %>%
#   mutate(Split = ifelse(any(atri_samples %in% semiVector(Split_1)),"Split_1",
#                         ifelse(any(atri_samples %in% semiVector(Split_2)),"Split_2",
#                                ifelse(any(atri_samples %in% semiVector(Split_3)),"Split_3","Split_4")))) %>%
#   ungroup() %>%
#   mutate(Species = "atri") %>%
#   select(Species,Loc,Split) %>%
#   write_tsv('data/output/Fixed_Locs/Atri.tsv')
# 
# acau_fixed_locs %>%
#   rowwise() %>%
#   mutate(Split = ifelse(any(acau_samples %in% semiVector(Split_1)),"Split_1",
#                         ifelse(any(acau_samples %in% semiVector(Split_2)),"Split_2",
#                                ifelse(any(acau_samples %in% semiVector(Split_3)),"Split_3","Split_4")))) %>%
#   ungroup() %>%
#   mutate(Species = "acau") %>%
#   select(Species,Loc,Split) %>%
#   write_tsv('data/output/Fixed_Locs/Acau.tsv')


# Read in all fixed sites
atub_fixed <- read_tsv('data/output/Fixed_Locs/Atub.tsv') %>% rename(Atub = "Split") %>% select(-Species)
awat_fixed <- read_tsv('data/output/Fixed_Locs/Awat.tsv') %>% rename(Awat = "Split") %>% select(-Species)
aret_fixed <- read_tsv('data/output/Fixed_Locs/Aret.tsv') %>% rename(Aret = "Split") %>% select(-Species)
apow_fixed <- read_tsv('data/output/Fixed_Locs/Apow.tsv') %>% rename(Apow = "Split") %>% select(-Species)
ahyp_fixed <- read_tsv('data/output/Fixed_Locs/Ahyp.tsv') %>% rename(Ahyp = "Split") %>% select(-Species)
aare_fixed <- read_tsv('data/output/Fixed_Locs/Aare.tsv') %>% rename(Aare = "Split") %>% select(-Species)
apal_fixed <- read_tsv('data/output/Fixed_Locs/Apal.tsv') %>% rename(Apal = "Split") %>% select(-Species)
atri_fixed <- read_tsv('data/output/Fixed_Locs/Atri.tsv') %>% rename(Atri = "Split") %>% select(-Species)
acau_fixed <- read_tsv('data/output/Fixed_Locs/Acau.tsv') %>% rename(Acau = "Split") %>% select(-Species)

# Get pure singletons

### Tub/Are ###

# 7,840
pure_atub_snps <- ind_signal %>%
  filter(Split_1 == atub_string | Split_2 == atub_string | Split_3 == atub_string | Split_4 == atub_string) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "Atub")

# 32,594
pure_aare_snps <- ind_signal %>%
  filter(Split_1 == aare_string | Split_2 == aare_string | Split_3 == aare_string | Split_4 == aare_string) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "Aare")

### Pow/Ret ###

# 3.266
pure_aret_snps <- ind_signal %>%
  filter(Split_1 == aret_string | Split_2 == aret_string | Split_3 == aret_string | Split_4 == aret_string) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "Aret")

# 13,233
pure_apow_snps <- ind_signal %>%
  filter(Split_1 == apow_string | Split_2 == apow_string | Split_3 == apow_string | Split_4 == apow_string) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "Apow")

### Hyp/Cau ###

# 10,091
pure_ahyp_snps <- ind_signal %>%
  filter(Split_1 == ahyp_string | Split_2 == ahyp_string | Split_3 == ahyp_string | Split_4 == ahyp_string) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "Ahyp")

# 19,413
pure_acau_snps <- ind_signal %>%
  filter(Split_1 == acau_string | Split_2 == acau_string | Split_3 == acau_string | Split_4 == acau_string) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "Acau")

### Pal/Wat/Spi ###

# 606 
pure_aspi_snps <- ind_signal %>% filter(str_detect(Singleton_Taxa,"trimmed_22002D-02-08_S35_L003_Nuclear")) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "Aspi")

# 979
pure_apal_snps <- ind_signal %>%
  filter(Split_1 == apal_string | Split_2 == apal_string | Split_3 == apal_string | Split_4 == apal_string) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "Apal")

# 2,378
pure_awat_snps <- ind_signal %>%
  filter(Split_1 == awat_string | Split_2 == awat_string | Split_3 == awat_string | Split_4 == awat_string) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "Awat")

# 10,450
palwat_string <- species_id_df %>% filter(Species_Pool %in% c("Apal",'Awat')) %>% pull(Sample_ID) %>% sort() %>% vectorSemi()
pure_palwat_snps <- ind_signal %>%
  filter(Split_1 == palwat_string | Split_2 == palwat_string | Split_3 == palwat_string | Split_4 == palwat_string) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "PalWat")

# 88,656
pure_atri_snps <- ind_signal %>%
  filter(Split_1 == atri_string | Split_2 == atri_string | Split_3 == atri_string | Split_4 == atri_string) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  mutate(Pure_Singleton = "Atri")


### Tub vs. Are ###
tub_are_splits <- full_join(atub_fixed,aare_fixed,by="Loc") %>% filter(!is.na(Atub) & !is.na(Aare)) %>% filter(Atub != Aare) %>% mutate(Split = "Atub_Aare") %>%
  filter(!Loc %in% pure_atub_snps$Loc & !Loc %in% pure_aare_snps$Loc)

tubare_split_contigs <- tub_are_splits %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  group_by(Contig) %>%
  summarize(TubAre_Split = n())

atub_diagnostic <- tubare_split_contigs %>%
  full_join(pure_atub_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_Atub = "n")) %>%
  mutate(TubAre_Split = ifelse(is.na(TubAre_Split),0,TubAre_Split)) %>%
  mutate(Pure_Atub = ifelse(is.na(Pure_Atub),0,Pure_Atub)) %>%
  filter(TubAre_Split > 0 | Pure_Atub > 0) %>%
  rowwise() %>%
  mutate(Atub_SNPs = sum(TubAre_Split,Pure_Atub)) %>%
  ungroup()

aare_diagnostic <- tubare_split_contigs %>%
  full_join(pure_aare_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_Aare = "n")) %>%
  mutate(TubAre_Split = ifelse(is.na(TubAre_Split),0,TubAre_Split)) %>%
  mutate(Pure_Aare = ifelse(is.na(Pure_Aare),0,Pure_Aare)) %>%
  filter(TubAre_Split > 0 | Pure_Aare > 0) %>%
  rowwise() %>%
  mutate(Aare_SNPs = sum(TubAre_Split,Pure_Aare)) %>%
  ungroup()

### Lower Clade A (Apal,Awat,Aspi) ###
aspi_id <- "trimmed_22002D-02-08_S35_L003_Nuclear"

aspi_fixed <- ind_signal %>%
  filter(!str_detect(Non_Base_Taxa, aspi_id)) %>%
  mutate(Aspi = case_when(
    str_detect(Singleton_Taxa, aspi_id) ~ "Singleton",
    str_detect(Split_1, aspi_id) ~ "Split_1",
    str_detect(Split_2, aspi_id) ~ "Split_2",
    str_detect(Split_3, aspi_id) ~ "Split_3",
    str_detect(Split_4, aspi_id) ~ "Split_3",
    TRUE ~ NA)) %>%
  select(Loc, Aspi)

all_lower_clade_a_fixed <- full_join(awat_fixed,apal_fixed) %>% full_join(aspi_fixed) %>% na.omit() %>%
  separate(Loc,c("Contig","Pos"),"/",FALSE)

aspi_diagnostic <- all_lower_clade_a_fixed %>%
  filter(!Loc %in% pure_aspi_snps$Loc) %>%
  filter(Aspi != Awat & Aspi != Apal ) %>%
  group_by(Contig) %>%
  summarize(Aspi_Split = n()) %>%
  full_join(pure_aspi_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_Aspi = "n")) %>%
  mutate(Aspi_Split = ifelse(is.na(Aspi_Split),0,Aspi_Split)) %>%
  mutate(Pure_Aspi = ifelse(is.na(Pure_Aspi),0,Pure_Aspi)) %>%
  filter(Aspi_Split > 0 | Pure_Aspi > 0) %>%
  rowwise() %>%
  mutate(Aspi_SNPs = sum(Aspi_Split,Pure_Aspi)) %>%
  ungroup()

awat_diagnostic <- all_lower_clade_a_fixed %>%
  filter(!Loc %in% pure_awat_snps$Loc) %>%
  filter(Awat != Aspi & Awat != Apal ) %>%
  group_by(Contig) %>%
  summarize(Awat_Split = n()) %>%
  full_join(pure_awat_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_Awat = "n")) %>%
  mutate(Awat_Split = ifelse(is.na(Awat_Split),0,Awat_Split)) %>%
  mutate(Pure_Awat = ifelse(is.na(Pure_Awat),0,Pure_Awat)) %>%
  filter(Awat_Split > 0 | Pure_Awat > 0) %>%
  rowwise() %>%
  mutate(Awat_SNPs = sum(Awat_Split,Pure_Awat)) %>%
  ungroup()

apal_diagnostic <- all_lower_clade_a_fixed %>%
  filter(!Loc %in% pure_apal_snps$Loc) %>%
  filter(Apal != Aspi & Apal != Awat ) %>%
  group_by(Contig) %>%
  summarize(Apal_Split = n()) %>%
  full_join(pure_apal_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_Apal = "n")) %>%
  mutate(Apal_Split = ifelse(is.na(Apal_Split),0,Apal_Split)) %>%
  mutate(Pure_Apal = ifelse(is.na(Pure_Apal),0,Pure_Apal)) %>%
  filter(Apal_Split > 0 | Pure_Apal > 0) %>%
  rowwise() %>%
  mutate(Apal_SNPs = sum(Apal_Split,Pure_Apal)) %>%
  ungroup()

palwat_diagnostic <- all_lower_clade_a_fixed %>%
  filter(Apal == Awat & Apal != Aspi ) %>%
  group_by(Contig) %>%
  summarize(PalWat_Split = n()) %>%
  full_join(pure_palwat_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_PalWat = "n")) %>%
  mutate(PalWat_Split = ifelse(is.na(PalWat_Split),0,PalWat_Split)) %>%
  mutate(Pure_PalWat = ifelse(is.na(Pure_PalWat),0,Pure_PalWat)) %>%
  filter(PalWat_Split > 0 | Pure_PalWat > 0) %>%
  rowwise() %>%
  mutate(PalWat_SNPs = sum(PalWat_Split,Pure_PalWat)) %>%
  ungroup()

### Lower Clade B (Acau, Ahyp, Aret, Apow) ###

all_lower_clade_b_fixed <- full_join(acau_fixed,ahyp_fixed) %>% full_join(aret_fixed) %>% full_join(apow_fixed) %>% na.omit() %>%
  separate(Loc,c("Contig","Pos"),"/",FALSE)

acau_diagnostic <- all_lower_clade_b_fixed %>%
  filter(!Loc %in% pure_acau_snps$Loc) %>%
  filter(Acau != Ahyp & Acau != Aret & Acau != Apow) %>%
  group_by(Contig) %>%
  summarize(Acau_Split = n()) %>%
  full_join(pure_acau_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_Acau = "n")) %>%
  mutate(Acau_Split = ifelse(is.na(Acau_Split),0,Acau_Split)) %>%
  mutate(Pure_Acau = ifelse(is.na(Pure_Acau),0,Pure_Acau)) %>%
  filter(Acau_Split > 0 | Pure_Acau > 0) %>%
  rowwise() %>%
  mutate(Acau_SNPs = sum(Acau_Split,Pure_Acau)) %>%
  ungroup()

ahyp_diagnostic <- all_lower_clade_b_fixed %>%
  filter(!Loc %in% pure_ahyp_snps$Loc) %>%
  filter(Ahyp != Acau & Ahyp != Aret & Ahyp != Apow) %>%
  group_by(Contig) %>%
  summarize(Ahyp_Split = n()) %>%
  full_join(pure_ahyp_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_Ahyp = "n")) %>%
  mutate(Ahyp_Split = ifelse(is.na(Ahyp_Split),0,Ahyp_Split)) %>%
  mutate(Pure_Ahyp = ifelse(is.na(Pure_Ahyp),0,Pure_Ahyp)) %>%
  filter(Ahyp_Split > 0 | Pure_Ahyp > 0) %>%
  rowwise() %>%
  mutate(Ahyp_SNPs = sum(Ahyp_Split,Pure_Ahyp)) %>%
  ungroup()

apow_diagnostic <- all_lower_clade_b_fixed %>%
  filter(!Loc %in% pure_apow_snps$Loc) %>%
  filter(Apow != Acau & Apow != Aret & Apow != Ahyp) %>%
  group_by(Contig) %>%
  summarize(Apow_Split = n()) %>%
  full_join(pure_apow_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_Apow = "n")) %>%
  mutate(Apow_Split = ifelse(is.na(Apow_Split),0,Apow_Split)) %>%
  mutate(Pure_Apow = ifelse(is.na(Pure_Apow),0,Pure_Apow)) %>%
  filter(Apow_Split > 0 | Pure_Apow > 0) %>%
  rowwise() %>%
  mutate(Apow_SNPs = sum(Apow_Split,Pure_Apow)) %>%
  ungroup()

aret_diagnostic <- all_lower_clade_b_fixed %>%
  filter(!Loc %in% pure_aret_snps$Loc) %>%
  filter(Aret != Acau & Aret != Aret & Aret != Ahyp) %>%
  group_by(Contig) %>%
  summarize(Aret_Split = n()) %>%
  full_join(pure_aret_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_Aret = "n")) %>%
  mutate(Aret_Split = ifelse(is.na(Aret_Split),0,Aret_Split)) %>%
  mutate(Pure_Aret = ifelse(is.na(Pure_Aret),0,Pure_Aret)) %>%
  filter(Aret_Split > 0 | Pure_Aret > 0) %>%
  rowwise() %>%
  mutate(Aret_SNPs = sum(Aret_Split,Pure_Aret)) %>%
  ungroup()

### Upper Clade (Atri) ###
other_upper_clade_samples <- species_id_df %>% filter(Upper_Clade == 1) %>% filter(is.na(Species_Pool)) %>% pull(Sample_ID) %>% sort()
upper_clade_samples <- species_id_df %>% filter(Upper_Clade == 1) %>% pull(Sample_ID) %>% sort()
single_atri <- species_id_df %>% filter(Species_Pool == "Atri") %>% head(1) %>% pull(Sample_ID)

upper_clade_fixed <- ind_signal %>%
  rowwise() %>%
  filter(!any(upper_clade_samples %in% vectorSemi(Non_Base_Taxa)))

atri_diagnostic <- upper_clade_fixed %>%
  filter(
    all(atri_samples %in% semiVector(Split_1)) & !any(other_upper_clade_samples %in% semiVector(Split_1)) |
      all(atri_samples %in% semiVector(Split_2)) & !any(other_upper_clade_samples %in% semiVector(Split_2)) |
      all(atri_samples %in% semiVector(Split_3)) & !any(other_upper_clade_samples %in% semiVector(Split_3)) |
      all(atri_samples %in% semiVector(Split_4)) & !any(other_upper_clade_samples %in% semiVector(Split_4))
  )

atri_splits <- atri_diagnostic %>%
  mutate(Atri = case_when(
    str_detect(Split_1,single_atri) ~ "Split_1",
    str_detect(Split_2,single_atri) ~ "Split_2",
    str_detect(Split_3,single_atri) ~ "Split_3",
    str_detect(Split_4,single_atri) ~ "Split_4",
    TRUE ~ NA    
  )) %>%
  select(Loc,Atri)

atri_diagnostic2 <- upper_clade_fixed %>%
  ungroup() %>%
  filter(!Loc %in% pure_atri_snps$Loc) %>%
  filter(Loc %in% atri_splits$Loc) %>%
  separate(Loc,c("Contig","Pos"),sep="/",FALSE) %>%
  group_by(Contig) %>%
  summarize(Atri_Split = n()) %>%
  full_join(pure_atri_snps %>% group_by(Contig) %>% count(Pure_Singleton) %>% ungroup() %>% select(Contig,Pure_Atri = "n")) %>%
  mutate(Atri_Split = ifelse(is.na(Atri_Split),0,Atri_Split)) %>%
  mutate(Pure_Atri = ifelse(is.na(Pure_Atri),0,Pure_Atri)) %>%
  filter(Atri_Split > 0 | Pure_Atri > 0) %>%
  rowwise() %>%
  mutate(Atri_SNPs = sum(Atri_Split,Pure_Atri)) %>%
  ungroup()

### All Diagnostic ###

major_locs <- bind_rows(upper_clade_locs,tubare_clade_locs,lower_clade_a_clade_locs,lower_clade_b_clade_locs,caud_hypo_clade_locs,pow_ret_clade_locs) %>%
  separate(Loc,c("Contig","Pos"),sep='/',remove = FALSE) %>%
  filter(Non_Base_Count < 5)

upper_clade_diagnostic <- major_locs %>%
  filter(SNP_Species == "Upper_Clade") %>%
  group_by(Contig) %>%
  summarize(Upper_Clade_SNP = n())

tubare_clade_diagnostic <- major_locs %>%
  filter(SNP_Species == "TubAre") %>%
  group_by(Contig) %>%
  summarize(TubAre_SNP = n())

lower_clade_a_diagnostic <- major_locs %>%
  filter(SNP_Species == "Lower_Clade_A") %>%
  group_by(Contig) %>%
  summarize(Lower_Clade_A_SNP = n())

lower_clade_b_diagnostic <- major_locs %>%
  filter(SNP_Species == "Lower_Clade_B") %>%
  group_by(Contig) %>%
  summarize(Lower_Clade_B_SNP = n())

diagnostic_list <- list(
  upper_clade_diagnostic,
  tubare_clade_diagnostic,
  lower_clade_a_diagnostic,
  lower_clade_b_diagnostic,
  atri_diagnostic2,
  atub_diagnostic,
  aare_diagnostic %>% select(-TubAre_Split),
  ahyp_diagnostic,
  acau_diagnostic,
  apow_diagnostic,
  aret_diagnostic,
  apal_diagnostic,
  awat_diagnostic,
  palwat_diagnostic,
  aspi_diagnostic
)

all_diagnostic <- reduce(diagnostic_list, full_join, by = "Contig") %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

sort_diagnostic <- all_diagnostic %>%
  select(Contig,Upper_Clade_SNP,TubAre_SNP,Lower_Clade_A_SNP,Lower_Clade_B_SNP,
         Atri_SNPs,
         Atub_SNPs,Aare_SNPs,
         Acau_SNPs,Ahyp_SNPs,Aret_SNPs,Apow_SNPs,
         Awat_SNPs,Apal_SNPs,PalWat_SNPs,Aspi_SNPs)


all_major_snps <- sort_diagnostic %>% 
  select(Contig, Upper_Clade_SNP, TubAre_SNP, Lower_Clade_A_SNP, Lower_Clade_B_SNP) %>%
  pivot_longer(
    cols = c(Upper_Clade_SNP, TubAre_SNP, Lower_Clade_A_SNP, Lower_Clade_B_SNP),names_to = "Clade",values_to = "Count") %>%
  group_by(Contig) %>%
  summarize(
    Total = sum(Count, na.rm = TRUE),
    Distinct = n_distinct(Clade[Count > 0 & !is.na(Count)])) %>%
  ungroup() %>%
  filter(Distinct == 4)

all_major_snp_contigs <- all_major_snps %>%
  pull(Contig)

tubpalwat_snps <- sort_diagnostic %>% 
  select(Contig,Lower_Clade_A_SNP, Apal_SNPs,Awat_SNPs,PalWat_SNPs,Atub_SNPs) %>%
  pivot_longer(
    cols = c(Lower_Clade_A_SNP, Apal_SNPs,Awat_SNPs,PalWat_SNPs,Atub_SNPs),names_to = "Clade",values_to = "Count") %>%
  group_by(Contig) %>%
  summarize(
    Total = sum(Count, na.rm = TRUE),
    Distinct = n_distinct(Clade[Count > 0 & !is.na(Count)])) %>%
  ungroup()

# Contigs of interest
# SISRS_contig-1825000053 (776 bp): 2+ for all 4 major groups, 2 Atub, 2 PalWat, 5 Atri, 1 Aspi
# SISRS_contig-632000056 (294 bp): 2 Lower Clade B (Pal/Wat/Spi), 2 Atub, 2 Awat, 1 Apal, 5 PalWat, 2 Acau
# SISRS_contig-917000029 (831 bp): 1 Lower Clade A (Cau/Hyp/Pow/Ret), 9 Atri, 1 Atub, 2 Ahyp, 3 Apow, 3 Aret, 1 Aspi