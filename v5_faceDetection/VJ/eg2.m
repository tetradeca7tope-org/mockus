
% Convert an OpenCV classifier XML file to a Matlab file 
ConvertHaarcasadeXMLOpenCV('HaarCascades/haarcascade_frontalface_alt.xml'); 
% Detect Faces 
Options.Resize=false; 
thresholds = 90*ones(22,1);
ObjectDetection('Images/1.jpg','HaarCascades/haarcascade_frontalface_alt.mat', ... 
  Options, thresholds);

