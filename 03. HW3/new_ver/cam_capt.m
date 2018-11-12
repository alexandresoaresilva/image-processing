function I = cam_capt
    % % webcam setup
    webcamlist_conn = webcamlist;
    webcam_connected = webcam(webcamlist_conn{1});
    webcam_connected.Resolution ='640x480';
    % webcam_connected.Resolution = '320x240';
    prv = preview(webcam_connected);
    prv.Parent.Visible = 'off';
%     prv.Parent.Parent.Parent.Position(1:2) = [0.6 .5];
    prv.Parent.Visible = 'on';
    prompt = 'Press enter to capture a frame';
    x = input(prompt);
    I = snapshot(webcam_connected);
    delete(prv);
end