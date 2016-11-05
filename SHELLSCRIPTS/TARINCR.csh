
#http://www.unixmen.com/performing-incremental-backups-using-tar/

#tar --listed-incremental=snapshot.file -cvzf backup.tar.gz /home/sandeepc/ > & ! log &
#tar --listed-incremental=snapshot.file -cvzf backup1.tar.gz /home/sandeepc/ > & ! log &

#This should be the next one...
#tar --listed-incremental=snapshot.file -cvzf backup2.tar.gz /home/sandeepc/ > & ! log &

#This should be the next one...
tar --listed-incremental=snapshot.file -cvzf backup3.tar.gz /home/sandeepc/ > & ! log &
