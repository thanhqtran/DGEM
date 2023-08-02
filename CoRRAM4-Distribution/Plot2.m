function Plot2(Data,nr,nc,fname,fsize,fmult,varargin)
%Plot the data in Data in nr by nc subplots
%{
    Author:     Alfred Maußner
    Purpose:    Plot the data in the structure Data in nr*nc
                subplots.

    Input:      Data, structure with length nr*nc, where 

                Data(i).l stores the lines to plot on the left y axis, 
                Data(i).r stores those to plot on the right y axis, and
                Data(i).n stores the legend names 
                nr, integer, the number of rows
                nc, integer, the number of columns
             fname, string, the name of the font used for axes numbers, axes labels, legends and titles
             fsize, the font size
             fmult, the multiplier of the font size used for axes labels (fmult(1)),
                    the legend (fmult(2)), and the title (fmult(3))

             Name, Value-Pairs
              'Box',  if present with value 1: show box around figure
             'Grid',  if present with value 1: show grid
            'Titles',  if present, cell of length nr*nc with optional titles for each subplot
          'AxLabels',  if present, cell of length 3, the x-axis, the left, and the right y-axis labels.

     Revision history:

            19      May 2019, first version, from Plot1.m
            20 November 2022, axis labels added for increased functionality

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

fh=figure;
fh.Units='centimeters';
fh.Position=[0.5 1 19 18];

[nobs, ~]=size(Data(1).l);

for i1=1:np
    if isempty(Data(i1).l)
        break;
    end
    subplot(nr,nc,i1);
    if ~isempty(Data(i1).r)
        yyaxis left;
        hl=plot(1:nobs,Data(i1).l);
        axl=gca;
        axl.FontName=fname;
        axl.FontSize=fsize;
        axl.LabelFontSizeMultiplier=fmult(1);        
        axl.YColor=myColor(1,:);
        if ~isempty(p.Results.AxLabels)
            xlabel(p.Results.AxLabels{1});
            ylabel(p.Results.AxLabels{2});
        else
            xlabel('Quarter');
            ylabel('Percent');        
        end        
        for j2=1:1:length(hl)
            hl(j2).LineStyle=myLineStyle{j2};
            hl(j2).LineWidth=1;
            hl(j2).Color=myColor(j2,:);
        end
        yyaxis right;
        hr=plot(1:nobs,Data(i1).r);
        axr=gca;
        axr.FontName=fname;
        axr.FontSize=fsize;
        axr.LabelFontSizeMultiplier=fmult(1);
        if ~isempty(p.Results.AxLabels)
            ylabel(p.Results.AxLables{3});
        else
            ylabel('Percent');
        end
        j3=length(hl);        
        axr.YColor=myColor(j3+1,:);
        for j2=1:1:length(hr)
            hr(j2).LineStyle=myLineStyle{j3+j2};
            hr(j2).LineWidth=1;
            hr(j2).Color=myColor(j3+j2,:);
        end
    else
        hl=plot(1:nobs,Data(i1).l);
        axl=gca;
        axl.FontName=fname;
        axl.FontSize=fsize;
        axl.LabelFontSizeMultiplier=fmult(1);
        if ~isempty(p.Results.AxLabels)
            xlabel(p.Results.AxLabels{1});
            ylabel(p.Results.AxLabels{2});
        else
            xlabel('Quarter');
            ylabel('Percent');        
        end        
        for j2=1:1:length(hl); hl(j2).LineStyle=myLineStyle{j2};hl(j2).LineWidth=1; hl(j2).Color=myColor(j2,:); end
    end
    axl.XLim=[1 nobs];
    if ~isempty(Data(i1).n)
        legend(Data(i1).n,'Location','northeast','FontName',fname,'FontSize',fsize*fmult(2));   
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
        title(p.Results.Titles{i1},'FontName',fname,'FontSize',fsize*fmult(3),'FontWeight','normal');
    end
end

return;   
end

