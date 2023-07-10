% 画出所选文件的STFT，并给出bin值
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\samples\';
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json', 2e6);
[loraSetSF6] = readLoraSet('sf6_BW31_25.json', 2e6);

% 生成preamble对应信道的idealchirp
% [downchirp, upchirp] = buildIdealchirp(loraSet, 375e3); % build idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0); % build idealchirp
% 读取文件夹下所有采样值文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));          
% for fileCount = 1:length(fileIn)
for fileCount = 1:1
    if mod(fileCount, 100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
    end
    [signal] = readSignalFile(inDir, fileIn(fileCount));
    % find Preamble
    [PreambleStartPos] = detect_preamble(signal, loraSet, downchirp, 1000); 
    if PreambleStartPos == 999
        continue;
    elseif PreambleStartPos ~= 1
        signal = circshift(signal, -(PreambleStartPos-1) * loraSet.dine);
    end
    % 调整cfo和timeoffset
    [cfo, windowsOffset] = get_cfo_winoff(signal, loraSet, downchirp, upchirp, loraSet.factor, true);
    signal = circshift(signal, -round(windowsOffset));
    [d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
%     [d_downchirp_cfoSF6, d_upchirp_cfoSF6] = rebuild_idealchirp_cfo(loraSetSF6, cfo, 0);
    stft(signal(1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     stft(signal(35.25*loraSet.dine:36.25*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
%     time_plot(signal(12.25*loraSet.dine+1:end), loraSet, d_upchirp_cfo, 5);
%     time_plot(signal(14.25*loraSet.dine+1:end), loraSet, d_downchirp_cfo, 80);
%     a = get_bin(signal(1:12*loraSet.dine), loraSet, d_downchirp_cfo, 10);
%     tmp = get_bin(signal(12.25*loraSet.dine+1 : end), loraSet, d_downchirp_cfo, 33);
%     disp((tmptmp - tmp)');
%     [singalOut] = divideChannel(loraSet, signal, [1,2,3,4], false);
%     loraSet.pass_arg = 0.05;
%     [singalOut] = lowPassFilterFir(signal, loraSet);
%     stft(singalOut(1, 1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
    writeSignalFile(loraSet, signal, 'E:\share\sf10\', 0, 0, 0, 10);
    
    
end

toc;
fclose all;