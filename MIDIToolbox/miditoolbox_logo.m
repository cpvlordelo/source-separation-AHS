%% MIDI toolbox 1.1 cover illustration% TE 18.5.2016a = readmidi('sarabande.mid');              % Sarabande for keysom, ivdistw = readmidi('wtcii01a.mid'); w=w(1:180,:); % wtcii01a for tempo curve%%clffigure(1)set(gcf,'Position', [440   190   706   608])FS=136;%%subplot(2,2,1) % som	keysom(a,0,'jet',14)	colormap jet	set(gca,'PlotBoxAspectRatio',[1 1 1])    ax=axis;    text(ax(2)/2+ax(2)*0.02,ax(4)/2+ax(4)*0.02,'m','FontSize',FS,'color','k','HorizontalAlignment','center','VerticalAlignment','middle')    text(ax(2)/2,ax(4)/2,'m','FontSize',FS,'color','w','HorizontalAlignment','center','VerticalAlignment','middle')%%subplot(2,2,2) % pcdist2	plotdist(ivdist2(a))    colorbar off	set(gca,'PlotBoxAspectRatio',[1 1 1])	colormap hot			colormap('hot(256)'); 		%map=colormap;		%newmap = flipud(map);		%colormap(newmap)         ax=axis;    text(ax(2)/2+ax(2)*0.02,ax(4)/2+ax(4)*0.02,'i','FontSize',FS,'color','k','HorizontalAlignment','center','VerticalAlignment','middle')    text(ax(2)/2,ax(4)/2,'i','FontSize',FS,'color','w','HorizontalAlignment','center','VerticalAlignment','middle')   %%subplot(2,2,3) % pianoroll	pianoroll(a(1:16,:),'hold','b')    set(gca,'Color',[.3 .7 .4])	set(gca,'PlotBoxAspectRatio',[1 1 1])         ax=axis;    text(ax(2)/2+ax(2)*0.01,ax(3)+(ax(4)-ax(3))/2+ax(4)*0.01,'d','FontSize',FS,'color','k','HorizontalAlignment','center','VerticalAlignment','middle')    text(ax(2)/2+ax(2)*0.02,ax(3)+(ax(4)-ax(3))/2+ax(4)*0.0,'d','FontSize',FS,'color','w','HorizontalAlignment','center','VerticalAlignment','middle')%%        subplot(2,2,4) % meter autocorr.    plot(onset(w,'sec'),tempocurve(w),'k','LineWidth',1.5)    xlabel('Time (s)'); ylabel('BPM')	set(gca,'PlotBoxAspectRatio',[1 1 1])    set(gca,'Color',[1 0 0])    ax=axis;    text(ax(2)/2+ax(2)*0.01,ax(3)+(ax(4)-ax(3))/2+ax(4)*0.01,'i','FontSize',FS,'color','k','HorizontalAlignment','center','VerticalAlignment','middle')    text(ax(2)/2+ax(2)*0.02,ax(3)+(ax(4)-ax(3))/2+ax(4)*0.0,'i','FontSize',FS,'color','w','HorizontalAlignment','center','VerticalAlignment','middle')%% write toolbox in the middletext(-45,92,'toolbox','FontSize',84,'FontName','Andale Mono')    %% printset(gcf,'PaperPositionMode','auto','InvertHardcopy','off','Color','w')print('miditoolbox_logo.png','-dpng','-r0')