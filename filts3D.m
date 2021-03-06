% HF should be your 4D AE signal
% ave and in should contain 3 values, the first value should be 0 or 1 to
% tell the function whether to try that type of filter, the 2nd should be
% the x direction value, the third should be z direction value

function X = filts3D(HF,ave,in,med,p,param)
%dims = [param.velmex.XNStep param.velmex.YNStep param.daq.HFdaq.pts param.daq.HFdaq.NoBurstTriggers];
dims = size(HF);
param = param;
if length(dims) ==3
    dims(4) = 1;
end

if med(1) == 1
    HF = real(HF);
    dims = size(HF);
    if strcmp(param.window,'beautify')
        if ndims(HF) > 2
            if med(2) > 1 && med(3) > 1 && med(4) == 1
                for i = 1:dims(3)
                    HF(:,:,i) = medfilt2(HF(:,:,i),[med(3) med(2)]);
                end
            end
        else
            if med(2) > 1 && med(3) > 1

                HF = medfilt2(HF,[med(3) med(2)]);

            elseif med(2) > 1 && med(3) == 1
                for i = 1:dims(1)
                    HF(i,:) = medfilt1(HF(i,:),med(2));
                end
            elseif med(3) > 1 && med(2) == 1
                for i = 1:dims(2)
                    HF(:,i) = medfilt1(HF(:,i),med(3));
                end
            end
        end
    else
        for i = 2:4
            if mod(med(i),2) ~= 1
                errordlg('Need window size for each dimension to be odd')
                return
            end
        end
        Sx = zeros(dims);
        if length(dims) < 4
            dims(4) = 1;
        end
        % H = ones(ave(2),ave(3))/(ave(2)*ave(3));

        for i = 1:dims(4)
            Sx(:,:,:,i) = medfilt3(HF(:,:,:,i),med(2:end));
            multiWaitbar('Median Filtering',i/dims(4));
        end
        HF = Sx;
        if p(1) > 1
            x = HF;
            if mod(p(1),2) == 0
                errordlg('window must be odd')
            else
                for i = 1:size(x,1)
                    for j = 1:size(x,2)
                        for k = 1:size(x,3)
                            x2(i,j,k,:) = medfilt1(squeeze(x(i,j,k,:)),p(1));
                        end
                    end
                    multiWaitbar('Median Filtering in Time',i/size(x,1));
                end
            end
            HF = x2;
            clear x2;
        end
    end
end

if ave(1) == 1
    switch param.window
        case 'Box'
            HF = real(HF);
            dims = size(HF);
            for i = 2:4
                if mod(ave(i),2) ~= 1
                    errordlg('Need window size for each dimension to be odd')
                    return
                end
            end
            Sx = zeros(dims);
            if length(dims) < 4
                dims(4) = 1;
            end
            % H = ones(ave(2),ave(3))/(ave(2)*ave(3));
            
            for i = 1:dims(4)
                Sx(:,:,:,i) = imboxfilt3(HF(:,:,:,i),ave(2:end));
                multiWaitbar('Smoothing Signal',i/dims(4));
            end
            clear HF
            HF = Sx;
            if p(2) > 1
                x = HF;
                if mod(p(2),2) == 0
                    errordlg('window must be odd')
                else
                    a = 1;
                    b = (1/p(2))*ones(1,p(2));
                    for i = 1:size(x,1)
                        for j = 1:size(x,2)
                            for k = 1:size(x,3)
                                x2(i,j,k,:) = filter(b,a,squeeze(x(i,j,k,:)));
                            end
                        end
                        multiWaitbar('Smoothing in Time',i/size(x,1));
                    end
                end
                HF = x2;
                clear x2;
            end
        case 'Tri'
            HF = mean_filter(HF,ave,p,'Tri');
        case '1D Gauss'
          HF = mean_filter(HF,ave,p,'1D Gauss');
        case 'beautify'
%             for i = 2:length(ave)
%                 if ave(i) > 1
%                     b(:,i-1) = gausswin(ave(i));
%                 end
%             end
%             b(:,4) = gausswin(p(2));
            Sx = zeros(dims);
            if ave(2) > 1 % X Dimension
                c = gausswin(ave(2));
              %  c = b(:,1); %X Dimension
                for i = 1:dims(1)
                    Sx(i,:) = filtfilt(c,1,squeeze(HF(i,:)));
                    multiWaitbar(' Smoothing in X',i/dims(2));
                end
                HF = Sx;
            end
            if ave(3) > 1 % Y Dimension
               c = gausswin(ave(3));
                for i = 1:dims(2)
                    Sx(:,i) = filtfilt(c,1,squeeze(HF(:,i)));
                    multiWaitbar('Smoothing in Y',i/dims(1));
                end
                HF = Sx;
            end
            case 'ND Gauss'
            HF = real(squeeze(HF));
            dims = size(HF);
            
            for i = 1:dims(end)
                if length(dims) == 4
                    Sx(:,:,:,i) = imgaussfilt3(HF(:,:,:,i),ave(2));
                elseif length(dims) == 3
                    Sx(:,:,i) = imgaussfilt(HF(:,:,i),ave(2));
                end
                multiWaitbar('Smoothing Signal',i/dims(4));
            end
            if ndims(Sx) == 3
                permute(Sx,[1 4 2 3]);
            end
            clear HF
            HF = Sx;
            if p(2) > 1
                x = HF;
                b = gausswin(p(2));
                for i = 1:size(x,1)
                    for j = 1:size(x,2)
                        for k = 1:size(x,3)
                            x2(i,j,k,:) = filter(b,1,squeeze(x(i,j,k,:)));
                        end
                    end
                    multiWaitbar('Smoothing in Time',i/size(x,1));
                end
                HF = x2;
                clear x2;
            end
            
        case 'Hamming'
            HF = mean_filter(HF,ave,p,'Hamming');
    end
    
    
    
end

if in(1) == 1
    if strcmp(param.window,'beautify')
        [x, y] = meshgrid(1:dims(1), 1:dims(2));
        [x1, y1] = meshgrid(1:1/in(2):dims(1),1:1/in(3):dims(2));
        X = interp2(x,y,HF,x1,y1);
    else
        for i = 1:3
            if size(HF,i) == 1
                P(i) = 1;
            else
                P(i) = 0;
            end
            q = sum(P);

        end
        if q == 0
            [x, y, z] = meshgrid(1:dims(2),1:dims(1),1:dims(3));
            [x2, y2, z2] = meshgrid(1:(1/in(3)):dims(2),1:(1/in(2)):dims(1),1:(1/in(4)):dims(3));
            %Ix = zeros(size(x2));
            for i = 1:dims(4)
                if in(5) == 0
                    Ix(:,:,:,i) = interp3(x,y,z,HF(:,:,:,i),x2,y2,z2);
                elseif in(5) == 1
                    Ix(:,:,:,i) = squarify(HF(:,:,:,i),'m');
                end
                multiWaitbar('Interpolating',i/dims(4))
            end
            X = Ix;

        end
        if q > 0
            %         if length(in) > 4
            in = in([1 2 4 5]);
            %             dims3 = dims(3);
            %             dims2 = dims(2);
            %         else
            %             dims3 = dims(2);
            %             dims2 = 1;
            %         end
            [x, z] = meshgrid(1:dims(1),1:dims(3));
            [x2,z2] = meshgrid(1:(1/in(2)):dims(1),1:(1/in(3)):dims(3));
            %         hf = ones(size(HF),'gpuArray');
            %         hf = hf.*HF;
            tic
            for i = 1:dims(4)
                for j = 1:dims(2)
                    if in(4) == 0
                        Ix(:,j,:,i) = interp2(x,z,squeeze(HF(:,j,:,i))',x2,z2);
                    elseif in(4) == 1
                        Ix(:,j,:,i) = permute(squarify(HF(:,j,:,i),'m'),[1 3 2]);
                    end
                    multiWaitbar('Interpolating',i/dims(4));
                end
            end
            toc
            X = permute(Ix,[3 2 1 4]);
        end
        if p(3) < 1
            Q = 'Decimating';
        elseif p(3) > 1
            Q = 'Interpolating';
        end
        clear x2
        if p(3) ~= 1
            x = X;

            K = size(x,3);
            J = size(x,2);
            for i = 1:size(x,1)
                for j = 1:J
                    for k = 1:K
                        x2(i,j,k,:) = interp1(linspace(0,1,size(x,4)),squeeze(x(i,j,k,:)),linspace(0,1,size(x,4)*p(3)));
                    end
                end
                multiWaitbar([Q ' in Time'],i/size(x,1));
            end
            X = x2;

        end
    end

    %       b1 = round(dims(3)/10);
    %     b2 = round(dims(3)*0.9);
    %     X = real(20*log10(real(X)./max(max(max(real(X(:,:,b1:b2,:)))))));
else
    X = HF;
    %     b1 = round(dims(3)/10);
    %     b2 = round(dims(3)*0.9);
    %     X = real(20*log10(real(X)./max(max(max(real(X(:,:,b1:b2,:)))))));

end

multiWaitbar('CLOSEALL');

function [HF] = mean_filter(HF,ave,p,win)
    HF = real(HF);
    dims = size(HF);
    HFmax = max(abs(HF(:)));
    
    switch win
        case 'Tri'
            for i = 2:length(ave)
                if ave(i) > 1
                    b(:,i-1) = triang(ave(i));
                end
            end
            b(:,4) = triang(p(2));
        case '1D Gauss'
            for i = 2:length(ave)
                if ave(i) > 1
                    b(:,i-1) = gausswin(ave(i));
                end
            end
            b(:,4) = gausswin(p(2));
        case 'Hamming'
            for i = 2:length(ave)
                if ave(i) > 1
                    b(:,i-1) = hamming(ave(i));
                end
            end
            b(:,4) = hamming(p(2));
    end
    Sx = zeros(dims);
    if ave(2) > 1 % X Dimension
        c = b(:,1); %X Dimension
        for i = 1:dims(end)
            for j = 1:dims(2)
                for k = 1:dims(3)
                    Sx(:,j,k,i) = filtfilt(c,1,squeeze(HF(:,j,k,i)));
                    %                     S = conv(c,squeeze(HF(:,j,k,i)));
                    %                     Sx(:,j,k,i) = interp1(linspace(0,1,length(S)),S,linspace(0,1,dims(1)));
                    %                     Sx(:,j,k,i) = circshift(Sx(:,j,k,i),floor(length(c)),1);
                end
            end
            multiWaitbar(' Smoothing in X',i/dims(4));
        end
        HF = Sx;
        clear Sx;
    end
    
    if ave(3) > 1
        c = b(:,2);
        for i = 1:dims(4)
            for j = 1:dims(1)
                for k = 1:dims(3)
                    Sx(j,:,k,i) = filtfilt(c,1,squeeze(HF(j,:,k,i)));
                end
            end
            multiWaitbar(' Smoothing in Y',i/dims(4));
        end
        HF = Sx;
        clear Sx;
    end
    if ave(4) > 1
        c = b(:,3);
        for i = 1:dims(4)
            for j = 1:dims(1)
                for k = 1:dims(2)
                    Sx(j,k,:,i) = filtfilt(c,1,squeeze(HF(j,k,:,i)));
                    %                       S = conv(c,squeeze(HF(j,k,:,i)));
                    %                     Sx(j,k,:,i) = interp1(linspace(0,1,length(S)),S,linspace(0,1,dims(1)));
                    
                end
            end
            multiWaitbar(' Smoothing in Z',i/dims(4));
        end
        HF = Sx;
        clear Sx
    end
    if p(2) > 1
        x = HF;
        c = b(:,4);
        for i = 1:size(x,1)
            for j = 1:size(x,2)
                for k = 1:size(x,3)
                    x2(i,j,k,:) = filter(c,1,squeeze(x(i,j,k,:)));
                end
            end
            multiWaitbar('Smoothing in Time',i/size(x,1));
        end
        HF = x2;
        clear x2;
    end
    HFmax2 = max(abs(HF(:)));
    HF = HF./(HFmax2/HFmax);




