function scroll(hTable, row)
    % Scroll the uitable to the specified row
    jScroll = findjobj(hTable);
    jScroll = jScroll.getViewport.getView;
    rowHeight = jScroll.getRowHeight;
    jScroll.scrollRectToVisible(java.awt.Rectangle(0, rowHeight * (row - 1), 1, rowHeight));
end
