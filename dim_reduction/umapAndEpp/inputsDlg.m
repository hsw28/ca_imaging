function [inputs, jd, jls, jts]=inputsDlg(msgTxt, title, labels, values, ...
    where, isNumeric, maxWidth, alignment, columns, ...
    low, high, nonModalCallback, javaWin, ...
    useApplyCloseBtns, mapFields, fieldCallback,...
    priorJd, sortLabels)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

if nargin<18
    sortLabels=true;
    if nargin<17
        priorJd=[];
        if nargin<16
            fieldCallback=[];
            if nargin<15
                mapFields=[];
                if nargin<14
                    useApplyCloseBtns=false;
                    if nargin<13
                        javaWin=[];
                        if nargin<12
                            nonModalCallback=[];
                            if nargin<11
                                high=[];
                                if nargin<10
                                    low=[];
                                    if nargin<9
                                        columns=2;
                                        if nargin<8
                                            %right=4, center=0, left=2
                                            alignment=2;
                                            if nargin<7
                                                maxWidth=14;
                                                if nargin<6
                                                    isNumeric=false;
                                                    if nargin<5
                                                        where='center';
                                                    end
                                                end
                                            end
                                            if isNumeric
                                                maxWidth=7;
                                                alignment=4;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
inputs={};
jd=[];
pnl=javaObjectEDT('javax.swing.JPanel', java.awt.BorderLayout(5,15));
if startsWith(msgTxt, '<html>')
    jl=javaObjectEDT('javax.swing.JLabel', msgTxt);
else
    app=BasicMap.Global;
    jl=javaObjectEDT('javax.swing.JLabel', ['<html>' app.h2Start ...
        msgTxt app.h2End '</html>']);
end
pnl.add(jl, 'North');
N=length(labels);

if sortLabels
    markerSorter=edu.stanford.facs.swing.MarkerSorter;
    I=markerSorter.sort1basedIndexes(labels);
    I(I==0)=length(labels);
    oLabels=labels(I);
    oValues=values(I);
else
    oValues=values;
    oLabels=labels;
    I=1:length(labels);
end
[south, jts, jls]=inputsPnl(oLabels,oValues,isNumeric,maxWidth,alignment,...
    java.awt.Color.BLUE, columns, [], [], low, high);
if ~isempty(mapFields)
    for j=1:N
        mapFields.set(Html.Remove(oLabels{j}), jts{j});
    end
end
if ~isempty(fieldCallback)
    for k=1:N 
        H=handle(jts{k}, 'CallbackProperties');
        set(H, ...
            'ActionPerformedCallback', @(h,e)doInput(h, k));
        set(H, 'FocusLostCallback', @(h,e)doInput(h, k));
    end
end
pnl.add(south, 'South');
msgType=javax.swing.JOptionPane.INFORMATION_MESSAGE;
pane=javaObjectEDT('javax.swing.JOptionPane', pnl, msgType);
pane.setOptionType(javax.swing.JOptionPane.OK_CANCEL_OPTION);
pane.setIcon(Gui.Icon('facs.gif'));
scs=get(0, 'ScreenSize');
d=south.getPreferredSize;
maxW=scs(3)*.8;
maxH=scs(4)*.5;
if d.width>maxW|| d.height>maxH
    dd=d;
    if d.width>maxW
        dd.width=maxW;
    end
    if d.height>maxH
        dd.height=maxH;
    end
    % add 25 pixesls for scrollbar
    dd.width=dd.width+25;
    south.setPreferredSize(dd);
end
isModal=isempty(nonModalCallback);
MatBasics.RunLater(@(h,e)focus(), .2);
if isempty(priorJd) || isModal
    jd=PopUp.Pane(pane, title, where, javaWin, isModal);
else
    jd=priorJd;
    jd.getContentPane.removeAll;
    jd.getContentPane.add(pane);
    jd.pack;    
end
if isModal
    paneValue=pane.getValue;
    if paneValue==0
        doInputs([], 'Ok');
    else
        doInputs([], 'Cancel');
    end
else
    Gui.SetDlgCloser(jd, @doInputs,{},useApplyCloseBtns)
    jd.setResizable(true);
end

    function focus
        jts{1}.requestFocus;
        jts{1}.selectAll
        drawnow;
    end

    function doInput(field, idx)
        if isequal(field.getForeground, Gui.ERROR_COLOR)
           return;
        end
        k=I(idx);
        old=values{k};
        if islogical(old)
            value=field.isSelected;
        else
            value=char(field.getText);
        end
        if ~isequal(value, old)
            values{k}=value;
            label=labels{k};
            fprintf('%s old=%s new=%d\n', Html.Remove(label), value, old);
            if isNumeric
                value=str2double(value);
            end
            feval(fieldCallback, value, label);
            jd.requestFocus;
        end
    end

    function ok=doInputs(~, btnLabel)
        ok=~strcmpi(btnLabel, 'Cancel');
        inputs={};
        if ok
            nV=length(oValues);
            inputs=cell(1, nV);
            for i=1:nV
                if islogical(oValues{i})
                    v=jts{i}.isSelected;
                else
                    v=char(jts{i}.getText);
                    badLimit=isequal(jts{i}.getForeground, Gui.ERROR_COLOR);
                    if badLimit
                        v=oValues{i};
                    end
                    if isNumeric
                        v=str2double(v);
                    end
                end
                k=I(i);
                inputs{k}=v;
            end
        end
        if ~isModal
            ok=feval(nonModalCallback, jd, btnLabel, inputs);
        end
    end
end