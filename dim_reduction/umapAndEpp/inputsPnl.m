function [pnl, jts, jls]=inputsPnl(labels,values, isNumeric, maxWidth, ...
    alignment, labelColor, columns, fncKeyTyped, fncFocusLost, ...
    low, high, tip)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

if nargin<11
    high=[];
    if nargin<10
        low=[];
        if nargin<9
            fncFocusLost=[];
            if nargin<8
                fncKeyTyped=[];
                if nargin<7
                    columns=1;
                    if nargin<6
                        labelColor=java.awt.Color.BLUE;
                        if nargin<5
                            %right=4, center=0, left=2
                            alignment=2;
                            if nargin<4
                                maxWidth=14;
                                if nargin<3
                                    isNumeric=false;
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
WEST_INSETS=java.awt.Insets(5,4,5,5);
EAST_INSETS=java.awt.Insets(5,5,5,10);

jls=cell(1, length(values));
jts=cell(1, length(values));
nV=length(values);
nL=length(labels);
if nL==0
    labels=cell(1,nV);
    for i=1:nV
        labels{i}=[num2str(i) '.'];
    end
    nL=nV;
end
pnl=javaObjectEDT('javax.swing.JPanel', java.awt.GridBagLayout);

c=java.awt.GridBagConstraints;
c.gridy=0;
c.gridx=0;
c.gridwidth=1;
gridx=0;
column=1;
for i=1:nV
    if i<=nL
        label=labels{i};
    else
        label=labels{end};
    end
    v=values{i};
    if isNumeric
        v=num2str(v);
    end
    c.gridx=gridx;
    if islogical(v)
        jt=javaObjectEDT('javax.swing.JCheckBox', label, v);
        c.gridwidth=2;
    else
        if ~isNumeric
            jl=javaObjectEDT('javax.swing.JLabel', label);
            jt=javaObjectEDT('javax.swing.JTextField', v, maxWidth);
            jt.setColumns(maxWidth);
            jt.setHorizontalAlignment(alignment);
            jj=handle(jt, 'CallbackProperties');
            if ~isempty(fncKeyTyped)
                set(jj, 'KeyTypedCallback', @(h,e)check(h,e,i));
            end
            if ~isempty(fncFocusLost)
                set(jj, 'FocusLostCallback', @(h,e)focusLost(h,e,i));
            end
            set(jj, 'FocusGainedCallback', @(h,e)focusGained(h,e));
        else
            [jt, jl]=Gui.AddNumberField(label, maxWidth, v, [], [], ...
                [], [], low, high);
        end
        jl.setForeground(labelColor);
        c.anchor=c.EAST;
        pnl.add(Gui.Panel(jl), c);
        jl.setHorizontalAlignment(javax.swing.JLabel.RIGHT);
        c.gridx=gridx+1;
        c.anchor=c.EAST;
    end
    jts{i}=jt;
    jls{i}=jl;
    pnl.add(Gui.Panel(jt),c);
    if column==columns
        c.gridy=c.gridy+1;
        gridx=0;
        column=1;
    else
        gridx=gridx+2;
        column=column+1;
    end
        
end
pnl=javaObjectEDT('javax.swing.JScrollPane', pnl);
    
    function check(h,e,idx)
        feval(fncKeyTyped, h, e, idx);
    end

    function focusLost (h, e, idx)
        feval(fncFocusLost, h, e, idx);
    end

    function focusGained(h, e)
        p=h.getParent;
        b=h.getBounds;
        p.scrollRectToVisible(b);
    end

end