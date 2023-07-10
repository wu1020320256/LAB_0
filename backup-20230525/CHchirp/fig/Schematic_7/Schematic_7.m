% 将输入文件夹中的信号文件读取分信道提取出来写入目的文件夹中
fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

% 采样值文件读取路径和保存路径
inDir = 'C:\Users\ZKevin\Desktop\CHchirp\fig\Schematic_7\';
% 设置preamble信道可选项
preambleChannelChoice = [-281.250e3, -93.750e3, 93.750e3, 281.250e3];
% 设置实验参数
PreambleChannel = 1;
payloadStartOffset = 13.25;
subchirpNum = 1;  % subchirp数目
payloadNum = 40;
channelNum = 4;
channelChoiceIndexTable = [1,2];
channelChoiceNum = length(channelChoiceIndexTable);  % 实际选择的信道数目
deleteFlag = false;
% 读取配置和验证文件
[loraSet] = readLoraSet('sf7_BW125.json');

% 生成0频段的idealchirp
[downchirp, upchirp] = buildIdealchirp(loraSet, 0);
% 读取文件夹下所有采样值,bin,channel文件
fileIn = dir(fullfile(inDir, '*.sigmf-data'));     
binTxt = dir(fullfile(inDir, 'bin_*.txt')); 
channelTxt = dir(fullfile(inDir, 'channel_*.txt'));
cfoTxt = dir(fullfile(inDir, 'cfo_*.txt'));
% 对所有文件按创建时间排序
[resultSort] = sortFileByTime([fileIn, binTxt, channelTxt, cfoTxt]);

dine = loraSet.dine;
fftX = loraSet.fft_x;
for fileCount = 1:length(fileIn)
    % 从文件中读取信号
    [signal] = readSignalFile(inDir, fileIn(resultSort(1, fileCount)));
    % 将信号按照信道存放规则划分
    signal = reshape(signal, [], 4).';
    % 调整信号的cfo和timeoffset
    cfo = load(strcat(inDir, cfoTxt(resultSort(4, fileCount)).name));

    % 根据cfo生成对应0频段的idealchirp
    [downchirpCfo, upchirpCfo] = rebuild_idealchirp_cfo(loraSet, cfo, 0);
    samples = signal(2, 13.25*dine+1:14.25*dine);
    dechirp = samples(1:dine/2) .* downchirpCfo(1:dine/2);
    dechirp_fft = abs(fft(dechirp, dine));
    dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
    figure(1);
    plot(dechirp_fft,'-','LineWidth',6); hold on;
    plot([64, 62], [2923, 259],'^','MarkerFaceColor','red','MarkerSize',20); hold on;
    xlim([0 fftX]);
    ylim([0 3000]);
    xlabel('BIN');
    ylabel('|FFT|(*10^3)'); 
    set(gca,'YTickLabel',[0 1 2 3]);
    set(gca,'FontSize',50,'FontName','Times New Roman');
    set(get(gca,'YLabel'),'FontSize',50);
    txt = {'BIN:64','|FFT|:2923'};
    text(47,2700,txt,'FontSize',50,'HorizontalAlignment','center','FontName','Times New Roman');
    txt = {'BIN:62','|FFT|:259'};
    text(43,500,txt,'FontSize',50,'HorizontalAlignment','center','FontName','Times New Roman');

%     dechirp = samples(1:dine/4) .* downchirpCfo(1:dine/4);
%     dechirp_fft = abs(fft(dechirp, dine));
%     dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
%     figure(2);
%     plot(dechirp_fft,'-','LineWidth',6); hold on;
%     plot([64, 60], [1648, 43],'^','MarkerFaceColor','red','MarkerSize',20); hold on;
%     xlim([0 fftX]);
%     ylim([0 3000]);
%     xlabel('BIN');
%     ylabel('|FFT|(*10^3)'); 
%     set(gca,'YTickLabel',[0 1 2 3]);
%     set(gca,'FontSize',50,'FontName','Times New Roman');
%     set(get(gca,'YLabel'),'FontSize',50);
%     txt = {'BIN:64','|FFT|:1648'};
%     text(47,1700,txt,'FontSize',50,'HorizontalAlignment','center','FontName','Times New Roman');
%     txt = {'BIN:60','|FFT|:43'};
%     text(43,500,txt,'FontSize',50,'HorizontalAlignment','center','FontName','Times New Roman');

%     dechirp = samples(1:dine/8) .* downchirpCfo(1:dine/8);
%     dechirp_fft = abs(fft(dechirp, dine));
%     dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
%     figure(3);
%     plot(dechirp_fft,'-','LineWidth',6); hold on;
%     plot([64, 56], [850, 15],'^','MarkerFaceColor','red','MarkerSize',20); hold on;
%     xlim([0 fftX]);
%     ylim([0 3000]);
%     xlabel('BIN');
%     ylabel('|FFT|(*10^3)'); 
%     set(gca,'YTickLabel',[0 1 2 3]);
%     set(gca,'FontSize',50,'FontName','Times New Roman');
%     set(get(gca,'YLabel'),'FontSize',50);
%     txt = {'BIN:64','|FFT|:850'};
%     text(64,1300,txt,'FontSize',50,'HorizontalAlignment','center','FontName','Times New Roman');
%     txt = {'BIN:56','|FFT|:15'};
%     text(43,500,txt,'FontSize',50,'HorizontalAlignment','center','FontName','Times New Roman');

%     dechirp = samples(1:dine/16) .* downchirpCfo(1:dine/16);
%     dechirp_fft = abs(fft(dechirp, dine));
%     dechirp_fft = dechirp_fft(1:fftX) + dechirp_fft(dine-fftX+1:dine);
%     figure(4);  
%     plot(dechirp_fft,'-','LineWidth',6); hold on;
%     plot([64, 48], [427, 13],'^','MarkerFaceColor','red','MarkerSize',20); hold on;
%     xlim([0 fftX]);
%     ylim([0 3000]);
%     xlabel('BIN');
%     ylabel('|FFT|(*10^3)'); 
%     set(gca,'YTickLabel',[0 1 2 3]);
%     set(gca,'FontSize',50,'FontName','Times New Roman');
%     set(get(gca,'YLabel'),'FontSize',50);
%     txt = {'BIN:64','|FFT|:427'};
%     text(64,800,txt,'FontSize',50,'HorizontalAlignment','center','FontName','Times New Roman');
%     txt = {'BIN:48','|FFT|:13'};
%     text(40,500,txt,'FontSize',50,'HorizontalAlignment','center','FontName','Times New Roman');

end
fclose all;