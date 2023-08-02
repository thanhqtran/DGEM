function Plot1(Data,nr,nc,fname,fsize,fmult,varargin)
%Plot the data in Data in nr by nc subplots
%{
    Author:     Alfred Maußner
    Purpose:    Plot the data in the structure Data in nr*nc
                subplots.

    Input: 
     Required:
        Data, structure with length nr*nc, where Data(i).x stores the lines to plot
              and Data(i).n stores the legends
          nr, integer, the number of rows
          nc, integer, the number of columns
       fname, string, the name of the font used for axes numbers, axes labels, legends and titles
       fsize, the font size
       fmult, the multiplier of the font size used for axes labels

    Name,Value-Pairs 
        'Box', if present with value 1: show box around figure
       'Grid', if present with value 1: show grid
     'Titles', if present, cell of length nr*nc with optional titles for each subplot
    'AxLabels', if present, cell of length 2, the x-axis and y-axsis labes

     Revision history:

            11   May 2019, first version
            28 March 2022, AxLabels added for increased functionality

 %}
p=inputParser;
p.addRequired('Data');
p.addRequired('nr');
p.addRequired('nc');
p.addRequired('fname');
p.addRequired('fsize');
p.addRequired('fmult');
p.addParameter('Box',true);
p.addParameter('Grid',true);
p.addParameter('Titles',[]);
p.addParameter('AxLabels',[]);
p.parse(Data,nr,nc,fname,fsize,fmult,varargin{:});

% Graphsettings
myColor=[[0 0 0];[0 0 1];[0 0.44 0];[1 0 0];[0 0.49 0.49];[0.46 0 0.46];[0.58 0.29 0];[0 0 0]];
myLineStyle={'-','-','-','-','-','-','-','--','--','--','--','--'};

np=length(Data);
if np>nr*nc
    error('There are more data than panels.');
    return; %#ok<UNRCH>
end
    
for i1=1:np
    if isempty(Data(i1).x)
        continue;
    end
    subplot(nr,nc,i1);
    [nobs, ~]=size(Data(i1).x);
    hl=plot(1:nobs,Data(i1).x);
    ax=gca;
    ax.FontName=fname;
    ax.FontSize=fsize;
    ax.LabelFontSizeMultiplier=fmult;
    ax.XLim=[1 nobs];
    if ~isempty(p.Results.AxLabels)
        xlabel(p.Results.AxLabels{1});
        ylabel(p.Results.AxLabels{2});
    else
        xlabel('Quarter');
        ylabel('Percent');        
    end
    for j2=1:1:length(hl); hl(j2).LineStyle=myLineStyle{j2};hl(j2).LineWidth=1; hl(j2).Color=myColor(j2,:); end
    if ~isempty(Data(i1).n)
        legend(Data(i1).n,'Location','northeast','FontName',fname,'FontSize',fsize*fmult);   
        if p.Results.Box
            legend('boxon');
        else
            legend('boxoff');
        end
    end
    if p.Results.Grid
        grid on;
    else
        grid off;
    end
    if ~isempty(p.Results.Titles)
        title(p.Results.Titles{i1},'FontName',fname,'FontSize',fsize*fmult,'FontWeight','normal');
    end
end

return;   
end

