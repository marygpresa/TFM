# TFM 

#configuro todo para que se automatice y suba a github (de esta manera no pierdo
# archivos y la gente puede ver, recomendar y replicar mi trabajo)
system("git config --global user.name 'marygpresa'")
setwd("/Users/mariagranados/TFM")  # Asegura que estás en el directorio correcto
system("git status")  # Comprueba si Git funciona aquí
pwd

commit_message <- paste("Actualización:", Sys.time())
system("git add .")  
system(paste("git commit -m '", commit_message, "'", sep=""))  
system("git push origin main")
