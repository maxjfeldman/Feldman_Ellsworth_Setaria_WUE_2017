

#####  START WITH TOP VIEW #####

##### tv_z3000 #####

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_tv_z3000_L2_MF.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_zoom_other \
-M imgtype:VIS,camera:TV,zoom:z3000 \
-T 40 \
-c

PID1=$!
wait $PID1

##### tv_z1000 #####

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_tv_z1000_L0_MF.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_zoom_other \
-M imgtype:VIS,camera:TV,zoom:z1000 \
-T 40

PID2=$!
wait $PID2

##### tv_z1_L2 #####

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_tv_z1_L2_MF.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-04-05-00-00_2014-02-08-03-00-00 \
-f imgtype_camera_zoom_other \
-M imgtype:VIS,camera:TV,zoom:z1 \
-T 40

PID3=$!
wait $PID3

##### tv_z1_L1 #####

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_tv_z1_L1_MF.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-08-05-00-00_2014-02-11-03-00-00-00 \
-f imgtype_camera_zoom_lifter_other \
-M imgtype:VIS,camera:TV,zoom:z1,lifter:h1 \
-T 40

PID4=$!
wait $PID4


##### tv_z1_L0 #####

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_tv_z1_L0_MF.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-11-03-00-00_2014-02-16-00-00-00-00 \
-f imgtype_camera_zoom_lifter_other \
-M imgtype:VIS,camera:TV,zoom:z1,lifter:h1 \
-T 40

PID5=$!
wait $PID5

#####  NOW SIDE VIEW #####

##### sv_z2500 #####
# Frame 0

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z2500_L2_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:0,zoom:z2500 \
-T 40

PID6=$!
wait $PID6


# Frame 90

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z2500_L2_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:90,zoom:z2500 \
-T 40

PID7=$!
wait $PID7


# Frame 180

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z2500_L2_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:180,zoom:z2500 \
-T 40

PID8=$!
wait $PID8

# Frame 270

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z2500_L2_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:270,zoom:z2500 \
-T 40

PID9=$!
wait $PID9


##### sv_z1000 #####
# Frame 0

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1000_L2_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:0,zoom:z1000 \
-T 40

PID10=$!
wait $PID10


# Frame 90

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1000_L2_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:90,zoom:z1000 \
-T 40

PID11=$!
wait $PID11


# Frame 180

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1000_L2_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:180,zoom:z1000 \
-T 40

PID12=$!
wait $PID12

# Frame 270

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1000_L2_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:270,zoom:z1000 \
-T 40

PID13=$!
wait $PID13

##### sv_z500 #####
# Frame 0

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z500_L1_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:0,zoom:z500 \
-T 40

PID14=$!
wait $PID14


# Frame 90

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z500_L1_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:90,zoom:z500 \
-T 40

PID15=$!
wait $PID15


# Frame 180

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z500_L1_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:180,zoom:z500 \
-T 40

PID16=$!
wait $PID16

# Frame 270

time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z500_L1_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:270,zoom:z500 \
-T 40

PID17=$!
wait $PID17

##### sv_z1_h1 #####

# Frame 0
time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1_L1_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-08-05-00-00-00_2014-02-11-03-00-00-00 \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:0,zoom:z1 \
-T 40

PID18=$!
wait $PID18

# Frame 90
time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1_L1_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-08-05-00-00-00_2014-02-11-03-00-00-00 \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:90,zoom:z1 \
-T 40

PID19=$!
wait $PID19

# Frame 180
time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1_L1_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-08-05-00-00-00_2014-02-11-03-00-00-00 \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:180,zoom:z1 \
-T 40

PID20=$!
wait $PID20

# Frame 270
time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1_L1_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-08-05-00-00-00_2014-02-11-03-00-00-00 \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:270,zoom:z1 \
-T 40

PID21=$!
wait $PID21


##### sv_z1_h0 #####

# Frame 0
time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1_L0_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-11-05-00-00-00_2014-02-16-00-00-00-00 \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:0,zoom:z1 \
-T 40

PID22=$!
wait $PID22

# Frame 90
time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1_L0_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-11-05-00-00-00_2014-02-16-00-00-00-00 \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:90,zoom:z1 \
-T 40

PID23=$!
wait $PID23

# Frame 180
time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1_L0_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-11-05-00-00-00_2014-02-16-00-00-00-00 \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:180,zoom:z1 \
-T 40

PID24=$!
wait $PID24

# Frame 270
time /home/mfeldman/plantcv_dev/plantcv-pipeline.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/tester/script/vis_sv_z1_L0_e82.py \
-s /home/mfeldman/tester/setaria_ril_db.sqlite3 \
-i /home/mfeldman/tester/out \
-D 2014-02-11-05-00-00-00_2014-02-16-00-00-00-00 \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:270,zoom:z1 \
-T 40

PID25=$!
wait $PID25
