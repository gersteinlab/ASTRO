# How to run the demo
## uncompress the gtf file
*tar -xJvf mmu.all.gtf.tar.xz*
## link a STAR reference
STAR reference is too big to be submitted to github. So you need downlaod or build a STAR reference of mouse (genome version:mm39). Then, change the starref variable in the json file to your own mm39 STAR index path.
## run the pipleine
Just run the following command in command line.  
*ASTRO test.json*
## Look at an example of output
Run 
*tar -xJvf output.tar.xz* 
will uncompress the output example.
But some buky files have been removed to reduce file size in the folder 
*rm output/STAR/tempLog.out; rm -rf output/STAR/temp_STARgenome/; rm -rf output/temps/barcode_db/*
