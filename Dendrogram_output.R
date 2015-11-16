
make.dendrogram(data = anim.data, regions = "AFR", variables = c("Mean.birth.rate", "Mean.vaccination"))

make.dendrogram(data = anim.data, regions = "AFR", variables = c("Coefficient.of.Variation",
                                                                 "Incidence"))

make.dendrogram(data = anim.data, regions = "AFR", variables = c("Coefficient.of.Variation",
                                                                 "Incidence",
                                                                 "Mean.birth.rate", 
                                                                 "Mean.vaccination"))



make.dendrogram(data = anim.data, regions = c("AMR","AFR"), variables = c("Mean.birth.rate", "Mean.vaccination"))

make.dendrogram(data = anim.data, regions = c("AMR","AFR"), variables = c("Coefficient.of.Variation",
                                                                          "Incidence"))

make.dendrogram(data = anim.data, regions = c("AMR","AFR"), variables = c("Coefficient.of.Variation",
                                                                          "Incidence",
                                                                          "Mean.birth.rate", 
                                                                          "Mean.vaccination"))
