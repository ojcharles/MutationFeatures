##### debug
# ensure files present in these directories
#ls /db
#ls /query
#ls /app



# Run R script
Rscript /tmp/.R 2>&1 >/dev/null

# cat the resultant file
cat /tmp/df.csv