% 将输入文件夹中的信号文件读取分信道提取出来写入目的文件夹中
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'E:\share\';
inDirVer = 'E:\uhd_share\';
% 设置preamble信道可选项
preambleChannelChoice = [-281.250e3, -93.750e3, 93.750e3, 281.250e3];
% 设置实验参数
preambleChannel = 1;
payloadStartOffset = 13.25;
subchirpNum = 4;
payloadNum = 40;
channelNum = 4;
channelChoiceIndexTable = [1,2,3,4];
channelChoiceNum = length(channelChoiceIndexTable);
% 读取配置和验证文件
[loraSet] = readLoraSet('sf10_BW125.json');

% 生成0频段的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0);
% 读取文件夹下所有采样值,bin,channel文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));     
binTxt = dir(fullfile(inDirVer, 'bin_*.txt')); 
channelTxt = dir(fullfile(inDirVer, 'channel_*.txt')); 
% 对所有文件按创建时间排序
[resultSort] = sortFileByTime([fileIn, binTxt, channelTxt]);

for fileCount = 1:length(fileIn)
    if mod(fileCount, 100) == 0           % 每循环100次，输出当前进度
        fprintf("The time is %d\n", fileCount);
    end
    % 从文件中读取信号文件
    [signal] = readSignalFile(inDir, fileIn(resultSort(1, fileCount)));
%     signal = reshape(signal, [], 4).';
%     signal = sum(signal, 1);
    % 调整信号的cfo和timeoffset
    [cfo, signal, packageFlag] = alignSignal(loraSet, signal, downchirp, upchirp, preambleChannelChoice(channelNum-preambleChannel+1));
    stft(signal(1:40*loraSet.dine), loraSet.sample_rate, 'Window',rectwin(64),'OverlapLength',32,'FFTLength',loraSet.fft_x);
    colormap(CustomColormap);
    colorbar('off');
    xlabel('Times(ms)');
    ylabel('Frequency(KHz)'); 
    set(gca,'FontSize',88);
    set(get(gca,'XLabel'),'FontSize',88);
    set(get(gca,'YLabel'),'FontSize',88);
    xlim([-inf 200]);
    ylim([-0.4 0.4001]);
    set(gca,'YTickLabel',[400 -200 0 200 400]);
    set(gca,'xtick',[],'xticklabel',[]);
    set(gca,'ytick',[],'yticklabel',[]);
    %     if packageFlag == false
%         continue;
%     end
%     % 将信号根据信道滤波划分
%     [singalOut] = divideChannel(loraSet, signal, preambleChannelChoice, false);
%     % 根据cfo生成对应0频段的idealchirp
%     [downchirpCfo, upchirpCfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
%     % 获得downchirp bin
%     [downchirpSync] = getDownchirpSync(loraSet, singalOut(preambleChannel,:), upchirpCfo);

end


fclose all;