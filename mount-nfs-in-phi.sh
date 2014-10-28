#!/bin/bash
mkdir -p /xeonphifs
chmod -R 777 /xeonphifs
/usr/sbin/exportfs -a -v
chkconfig nfs on
service nfs restart
scp root@mic0:/etc/fstab xeon-phi-fstab
echo "host:/xeonphifs /xeonphifs nfs rsize=8192,wsize=8192,nolock,intr 0 0" >> xeon-phi-fstab
scp xeon-phi-fstab root@mic0:/etc/fstab
ssh root@mic0 "mkdir /xeonphifs && mount -a"
