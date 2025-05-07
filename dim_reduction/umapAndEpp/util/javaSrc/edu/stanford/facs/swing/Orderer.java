package edu.stanford.facs.swing;

import java.awt.*;
import java.awt.dnd.*;
import java.awt.event.*;
import java.util.*;
import java.util.List;

import javax.swing.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.IOException;

import javax.swing.event.*;



public  class Orderer<T> {
  	
	private Window owner;
	final JCheckBox cbAutoUpdate=new JCheckBox("Update table now");
	private boolean showUpdate = true;
	private boolean emptyNotAllowed=false;
	public Orderer(final AbstractButton updateButton, final boolean autoUpdate, boolean showUpdate, boolean treeOnlyContext, boolean emptyNotAllowed){
		this(updateButton, autoUpdate, showUpdate, treeOnlyContext);
		this.emptyNotAllowed = emptyNotAllowed;
	}
	public Orderer(final AbstractButton updateButton, final boolean autoUpdate, boolean showUpdate, boolean treeOnlyContext){
		this(updateButton, autoUpdate, showUpdate);
		if (treeOnlyContext) {
			cbAutoUpdate.setText("Update tree now");
		}
	}
	public Orderer(final AbstractButton updateButton, final boolean autoUpdate, boolean showUpdate){
		this(updateButton, autoUpdate);
		this.showUpdate = showUpdate;
	}
	public Orderer(final AbstractButton updateButton, final boolean autoUpdate){
		this.updateButton=updateButton;
		cbAutoUpdate.setSelected(autoUpdate);		
  		if (updateButton != null){
  			updateButton.setEnabled(false);
  			updateButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					updateButton.setEnabled(false);
				}
			});
			updateButton.setEnabled(false);
			cbAutoUpdate.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					if (cbAutoUpdate.isSelected()){
						if (Orderer.this.updateButton.isEnabled()){
							updateButton.doClick();
						}
					}
				}
			});			
  		}
	}

	public void setButtonEnabled(boolean enabled){
		if (updateButton!=null){
			updateButton.setEnabled(enabled);
		}
	}
	private void enableUpdates(){
		if (updateButton!=null){
			updateButton.setEnabled(true);
			if(cbAutoUpdate.isSelected()){
				updateButton.doClick();
				SwingUtilities.invokeLater(new Runnable() {
					
					@Override
					public void run() {
						// TODO Auto-generated method stub
						updateButton.setEnabled(false);
						
					}
				});
			}
			ToolTipOnDemand.getSingleton().showLater(updateButton);
		}
	}

	public JDialog dlg;
	
	public JList<T> leftList, rightList;
	public JScrollPane leftListScrollPane, rightListScrollPane;
	private final AbstractButton updateButton;
	private DefaultListModel<T> leftModel, rightModel=new DefaultListModel<T>();
	
	public JLabel setCountListeningLabel(final JLabel lbl, final String prefixNoHtml, final String startHtml, final String endHtml){
		leftModel.addListDataListener(new ListDataListener() {
	
			public void intervalRemoved(ListDataEvent e) {
				System.out.println("intervalRemoved");
				reflectCount(lbl, prefixNoHtml, startHtml, endHtml);
			}
			
			public void intervalAdded(ListDataEvent e) {
				System.out.println("intervalAdded");
				reflectCount(lbl, prefixNoHtml, startHtml, endHtml);
			}
			
			public void contentsChanged(ListDataEvent e) {
				System.out.println("intervalRemoved");
				reflectCount(lbl, prefixNoHtml, startHtml, endHtml);
			}
		});
		return lbl;
	}
	public void reflectCount(final JLabel lbl, final String prefixNoHtml, final String startHtml, final String endHtml){
		final int sz=leftList.getModel().getSize();
		final int cnt=sz+rightList.getModel().getSize();
		lbl.setText("<html>"+prefixNoHtml+startHtml+sz+"/"+cnt+endHtml+"</html>");
	}
	public void resetRight(final Collection<T> items){
		final Set<T> set = new LinkedHashSet<T>(items);
		rightModel.removeAllElements();
		final Iterator<T> it=set.iterator();
		while(it.hasNext()){
			final T o=it.next();
			if (!leftModel.contains(o)){
				rightModel.addElement(o);
			}
		}
	}
	
	public void resetLeft(final Collection<T> items){
		final Set<T> set = new LinkedHashSet<T>(items);
		leftModel.removeAllElements();
		final Iterator<T> it=set.iterator();
		while(it.hasNext()){
			final T o=it.next();
			if (!rightModel.contains(o)){
				leftModel.addElement(o);
			}
		}
	}
	
	public static interface DoubleClickResponder{
		void respond(final Object clickedOn, final Component parent);
	}

	public static abstract class StringTransferHandler extends TransferHandler {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		protected abstract String exportString(JComponent c);
		protected abstract void importString(JComponent c, String str);
		protected abstract void cleanup(JComponent c, boolean remove);

		protected Transferable createTransferable(JComponent c) {
			return new StringSelection(exportString(c));
		}

		public int getSourceActions(JComponent c) {
			return MOVE;
		}

		public boolean importData(final JComponent c, final Transferable t) {
			if (canImport(c, t.getTransferDataFlavors())) {
				try {
					String str = (String)t.getTransferData(DataFlavor.stringFlavor);
					importString(c, str);
					return true;
				} catch (UnsupportedFlavorException ufe) {
				} catch (IOException ioe) {
				}
			}

			return false;
		}

		protected void exportDone(JComponent c, Transferable data, int action) {
			cleanup(c, action == MOVE);
		}

		public boolean canImport(JComponent c, DataFlavor[] flavors) {
			for (int i = 0; i < flavors.length; i++) {
				if (DataFlavor.stringFlavor.equals(flavors[i])) {
					return true;
				}
			}
			return false;
		}
	}
	
	public JLabel rightListHeader = new JLabel("<html><b>Available columns</html>");
	public void setRightListHeader(final String txt){
		rightListHeader.setHorizontalAlignment(JLabel.CENTER);
		rightListHeader.setText(txt);
	}
	
	public void resetLists(
			final Collection<T> leftItems,
			final Collection<T> rightItems) {
		leftModel.removeAllElements();
		Iterator<T> it=leftItems.iterator();
		while(it.hasNext()){
			leftModel.addElement(it.next());
		}
		rightModel.removeAllElements();
		it=rightItems.iterator();
		while(it.hasNext()){
			rightModel.addElement(it.next());
		}
	}
	
	public JPanel setLists(
			final Window owner,
			final Collection<T> leftItems,
			final Collection<T> rightItems,
			final JComponent atBottom,
			final JComponent leftLabel,
			final JComponent rightLabel,
			final ListCellRenderer<T> leftLcr,
			final ListCellRenderer<T> rightLcr,
			final DoubleClickResponder leftDoubleClick
			) {
			this.dlg = new JDialog(owner);
			this.owner=owner;
			dlg.setModal(true);
			final JPanel mainPanel=this.getMainPanel(dlg, 
					leftItems, rightItems, 
					atBottom, leftLabel, rightLabel, 
					leftLcr, rightLcr, leftDoubleClick);
			dlg.setContentPane(mainPanel);
			dlg.setTitle(title);
			dlg.pack();
			return mainPanel;
	}
	
	public JPanel getMainPanel(
			final Collection<T> leftList,
			final String leftListLabel,
			final Collection<T> rightList,
			final String rightListLabel){
		final JLabel leftLabel=new JLabel(leftListLabel), 
				rightLabel=new JLabel(rightListLabel);
		rightLabel.setHorizontalAlignment(JLabel.CENTER);
		leftLabel.setHorizontalAlignment(JLabel.CENTER);
		return this.getMainPanel(null, 
				leftList, rightList, null, 
				leftLabel, rightLabel, null,null,null);

	}
	
	public int visibleRows=0;
	public String itemName="keyword", itemNamePlural="keyword(s)";
	public boolean resetRequired= true;
	
	private void sort(DefaultListModel dlm) { 
	    Object [] dlma = dlm.toArray();   
	    List dlml=Arrays.asList(dlma);
	    if (ms==null) {
	    	Collections.sort(dlml, new Comparator<String>() {
	    		@Override
	    		public int compare(String o1, String o2) {
	    			if (o1.startsWith("$")){
	    				o1=o1.substring(1);
	    			}
	    			if (o2.startsWith("$")){
	    				o2=o2.substring(1);
	    			}
	    			return o1.compareTo(o2);
	    		}
	    	});
	    } else {
	    	dlml=ms.sort(dlml);
	    }
	    dlm.clear(); 
	    for (Object x : dlml) {
	    	dlm.addElement(x); 
	    }
	}
	
	private String rightListTip, leftListTip;
	private JButton leftButton, rightButton;
	private boolean isLeftRightOk=true;	
	public boolean setIsLeftRightOk(final boolean ok) {
		boolean prev=ok;
		this.isLeftRightOk=ok;
		if (this.leftButton != null) {
			this.leftButton.setVisible(this.isLeftRightOk);
		}
		if (this.rightButton != null) {
			this.rightButton.setVisible(this.isLeftRightOk);
		}
		if (this.leftList != null) {
			if (ok){
				leftList.setToolTipText(
						"<html>To <b><i>stop using</i></b> a "+itemName+" <u>drag</u> it right to the "+rightListTip+
						" list or click right button.<br>To <b><i>reorder</i></b> a "+itemName+" drag up/down or click up/down buttons.</html>");
			}else {
				leftList.setToolTipText(
						"<html>To <b><i>reorder</i></b> a "+itemName+" drag up/down or click up/down buttons.</html>");
				
			}
		}
		if (this.rightList != null) {
			if (ok) {
				rightList.setToolTipText(
					"<html>To <b><i>use</i></b> a "+itemName+" <u>drag</u> it left to the "+leftListTip+
					" list.</html>"); 
			}else {
				rightList.setToolTipText(
						"<html>To <b><i>reorder</i></b> a "+itemName+" drag up/down or click up/down buttons.</html>");
			}
		}
		return prev;
	}
	
	private JComponent lastDrag;
	private boolean dropOk=true;
	public JPanel getMainPanel(
			final Window wnd,
			final Collection<T> leftItems,
			final Collection<T> rightItems,
			final JComponent atBottom,
			final JComponent leftLabel,
			final JComponent rightLabel,
			final ListCellRenderer<T> leftLcr,
			final ListCellRenderer<T> rightLcr,
			final DoubleClickResponder doubleClickLeftList
			) {
		final JPanel buttonPanel=new JPanel();
		class ListItems {
			private Set<T> set=new TreeSet<T>();
			public Object stringToObject(String str){
				final Iterator<T> it=set.iterator();
				while(it.hasNext()){
					Object obj=it.next();
					if(obj.toString().equals(str)){
						return obj;
					}
				}
				return null;
			}
			public void addAll(Collection<T> data){
				set.addAll(data);
			}
		}
		final Set<T> leftSet = new LinkedHashSet<T> (leftItems);
		final Set<T> rightSet = new LinkedHashSet<T> (rightItems);
		final int rows;
		if (visibleRows==0){
			int n1=leftSet.size(),n2=rightSet.size();
			int n=n1+n2;
			if (n<5){
				n=5;
			} else if (n>15){
				n=15;
			}
			rows=n;
		}else{
			rows=visibleRows;
		}
		leftModel=new DefaultListModel<T>();
		final JDialog dlg=this.dlg;
		final FocusListener fl=new FocusListener() {
			public void focusLost(FocusEvent e) {
			}
			@Override
			public void focusGained(FocusEvent e) {
				if (dlg != null){
					dlg.getRootPane().setDefaultButton(null);
				}
			}
		};
		leftList = new JList<T>(leftModel){
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			public void setListData(final T[] listData) {
				DefaultListModel<T> model=new DefaultListModel<T>();
				for(int index=0;index<listData.length;index++){
					model.addElement(listData[index]);
				}
				setModel (model);
			}
		};
		if (leftLcr != null){
			leftList.setCellRenderer(leftLcr);
		}
		leftList.setBorder(BorderFactory.createEmptyBorder(10, 5, 8,5));
		leftList.addFocusListener(fl);

		final ListItems listItems=new ListItems();
		listItems.addAll(leftSet);
		listItems.addAll(rightSet);
		rightModel=new DefaultListModel<T>();
		Iterator<T> it=leftSet.iterator();
		while(it.hasNext()){
			leftModel.addElement(it.next());
		}
		it=rightSet.iterator();
		while(it.hasNext()){
			rightModel.addElement(it.next());
		}
		rightList = new JList<T>(rightModel){
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			public void setListData(final T[] listData) {
				DefaultListModel<T> model=new DefaultListModel<T>();
				for(int index=0;index<listData.length;index++){
					model.addElement(listData[index]);
				}
				setModel (model);
			}
		};
		
		rightList.addFocusListener(fl);
		if (rightLcr != null){
			rightList.setCellRenderer(rightLcr);
		}
		rightList.setBorder(BorderFactory.createEmptyBorder(10, 5, 8, 5));
		rightList.setVisibleRowCount(rows);
		rightList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		rightList.setDragEnabled(true);
		JPanel mainPanel = null;
		final DragSource ds=new DragSource();
		class ListDragGesture implements DragGestureListener,DragSourceListener{
			JComponent m_cmp=null;
			public ListDragGesture(final JComponent cmp){
				m_cmp=cmp;
				ds.createDefaultDragGestureRecognizer(cmp,
						DnDConstants.ACTION_COPY, this);
			}
			public void dragDropEnd(DragSourceDropEvent dsde) {}
			public void dragEnter(DragSourceDragEvent dsde) {}
			public void dragExit(DragSourceEvent dse) {}
			public void dragOver(DragSourceDragEvent dsde) {}
			public void dropActionChanged(DragSourceDragEvent dsde) {}
			public void dragGestureRecognized(final DragGestureEvent dge) {
				try{
					ds.startDrag(dge, DragSource.DefaultMoveDrop,				
							((StringTransferHandler) ((JList<?>)m_cmp).getTransferHandler())
							.createTransferable(m_cmp), this);
					lastDrag=m_cmp;					
				}
				catch(final Exception e){						
				}
			}
		}
		
		class ListTransferHandler extends StringTransferHandler {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;
			private int[] indices = null;
			private int addIndex = -1; // Location where items were added
			private int addCount = 0; // Number of items added.
			public ListTransferHandler(){
			}
			// Bundle up the selected items in the list
			// as a single string, for export.
			protected String exportString(final JComponent c) {
				JList<?> source = (JList<?>) c;
				indices = source.getSelectedIndices();
				@SuppressWarnings("deprecation")
				Object[] values = source.getSelectedValues();
				final StringBuilder buff = new StringBuilder();
				for (int i = 0; i < values.length; i++) {
					Object val = values[i];
					buff.append(val == null ? "" : val.toString());
					if (i != values.length - 1) {
						buff.append("\n");
					}
				}
				return buff.toString();
			}

			// Take the incoming string and wherever there is a
			// newline, break it into a separate item in the list.
			protected void importString(final JComponent c, final String str) {
				if(!(c instanceof JList) || str.trim().length()==0){
					return;
				}
				final boolean same=c.equals(lastDrag);
				dropOk=same || isLeftRightOk;
				System.out.println("importString() drag & drop are same = "+same + ", drop OK="+dropOk);
				if (!dropOk) {
					return;
				}
				@SuppressWarnings("unchecked")
				JList<T> target = (JList<T>) c;
				DefaultListModel<T> listModel=(DefaultListModel<T>)target.getModel();
				int index = target.getSelectedIndex()==-1? 0 : target.getSelectedIndex();
				addIndex = index;
				String[] values = str.split("\n");
				addCount = values.length;
				for (int i = 0; i < values.length; i++) {
					@SuppressWarnings("unchecked")
					T obj=(T)listItems.stringToObject(values[i]);
					if(rightList.equals(target)){
						((DefaultListModel<T>)leftList.getModel()).removeElement(obj);
					}
					if(leftList.equals(target)){
						((DefaultListModel<T>)rightList.getModel()).removeElement(obj);
					}
					listModel.removeElement(obj);
					try{
						listModel.add(index++,obj);
					}catch(Exception e){
						listModel.add(listModel.getSize()-1,obj);
					}
				}
				enableUpdates();
			}

			protected void cleanup(final JComponent c, final boolean remove) {
				final boolean same=c.equals(lastDrag);
				System.out.println("cleanup() drag & drop are same = "+same + ", drop OK="+dropOk);
				if (!dropOk) {
					return;
				}
				if (remove && indices != null) {
					if (addCount > 0) {
						for (int i = 0; i < indices.length; i++) {
							if (indices[i] > addIndex) {
								indices[i] += addCount;
							}
						}
					}
				}
				indices = null;
				addCount = 0;
				addIndex = -1;
			}
		}


		leftList.setTransferHandler(new ListTransferHandler());
		ds.createDefaultDragGestureRecognizer(leftList,
				DnDConstants.ACTION_MOVE, new ListDragGesture(leftList));
		rightList.setTransferHandler(new ListTransferHandler());
		ds.createDefaultDragGestureRecognizer(rightList,
				DnDConstants.ACTION_MOVE, new ListDragGesture(rightList));
		leftListTip="<b>"+(leftLabel instanceof JLabel ? ( (JLabel)leftLabel).getText() : " left list " )+"</b>";
		rightListTip="<b>"+(rightLabel instanceof JLabel ? ( (JLabel)rightLabel).getText() : " right list " )+"</b>";
		final JPanel moveButtonPanel2 = new JPanel(new BorderLayout());
		final JButton top2 = new ImageButton(SwingUtil.getImageGifIcon("moveTop"));
		ActionListener anAction = new ActionListener() {
			public void actionPerformed(final ActionEvent a) {
				Mover.move(rightList, Mover.Direction.TOP);
				rightList.requestFocus();
				enableUpdates();
			}
		};
		SwingUtil.echoAction(rightList, top2, anAction, 
				KeyStroke.getKeyStroke(KeyEvent.VK_HOME, InputEvent.ALT_MASK), 't');
		final JButton bottom2=new ImageButton(SwingUtil.getImageGifIcon("moveBottom"));
		anAction = new ActionListener() {
			public void actionPerformed(final ActionEvent a) {
				Mover.move(rightList, Mover.Direction.BOTTOM);
				rightList.requestFocus();
				enableUpdates();
			}
		};
		SwingUtil.echoAction(leftList, bottom2, anAction, 
				KeyStroke.getKeyStroke(KeyEvent.VK_END, InputEvent.ALT_MASK), 'b');
		
		final JButton moveUp2 = new JButton(SwingUtil.getImageGifIcon("upArrow"));
		moveUp2.setToolTipText("<html>Click to move a used "+itemName+" up the " + leftListTip + " list</html>");
		final JButton moveDown2=new JButton(SwingUtil.getImageGifIcon("downArrow"));
		moveDown2.setToolTipText("<html>Click to move a used "+itemName+" down the " + leftListTip + " list</html>");
		moveUp2.setMargin(new Insets(1, 2, 1, 2));
		moveDown2.setMargin(new Insets(1, 2, 1, 2));
		anAction = new ActionListener() {
			public void actionPerformed(final ActionEvent a) {
				Mover.move(rightList, Mover.Direction.UP);
				rightList.requestFocus();
				enableUpdates();
			}
		};
		SwingUtil.echoAction(rightList, moveUp2, anAction, 
				KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.ALT_MASK), 'u');
		anAction = new ActionListener() {
			public void actionPerformed(final ActionEvent a) {
				Mover.move(rightList, Mover.Direction.DOWN);
				rightList.requestFocus();
				enableUpdates();
			}
		};
		SwingUtil.echoAction(rightList, moveDown2, anAction, 
				KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.ALT_MASK), 'd');
		final JPanel northPanel2 = new JPanel(new BorderLayout(0, 4));
		northPanel2.add(top2, "North");
		northPanel2.add(moveUp2, "South");
		moveButtonPanel2.add(northPanel2, "North");
		final JPanel southPanel2 = new JPanel(new BorderLayout(0, 4));
		southPanel2.add(moveDown2, "North");
		southPanel2.add(bottom2, "South");
		moveButtonPanel2.add(southPanel2, "South");
		final ListSelectionListener selectAndOrderListSelectionListener2 = new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent evt) {
				if (rightList.getSelectedIndex() > -1) {
					leftButton.setEnabled(true);
					ToolTipOnDemand.getSingleton().showLater(rightButton, false, null, -35, 45, false, 
							"<html><center><u>Double click</u> to <b><i>stop using</i></b><br><small>(or <u>drag</u> to "+rightListTip+" list<br>or <u>click</u> right button)</small></center><</html>");
				} else {
					leftButton.setEnabled(false);
				}
				top2.setEnabled(true);
				bottom2.setEnabled(true);
				moveUp2.setEnabled(true);
				moveDown2.setEnabled(true);
				int[] indices=rightList.getSelectedIndices();
				for(int i=0;i<indices.length;i++){
					if(indices[i]==0){
						top2.setEnabled(false);
						moveUp2.setEnabled(false);
					}
					else if(indices[i]==rightList.getModel().getSize()-1){
						bottom2.setEnabled(false);
						moveDown2.setEnabled(false);
					}
				}
				if(rightList.getModel().getSize() <=1 || indices.length <= 0) {
					moveUp2.setEnabled(false);
					moveDown2.setEnabled(false);
					top2.setEnabled(false);
					bottom2.setEnabled(false);
				}
				else if(rightList.getSelectedIndices().length == rightList.getModel().getSize()) {
					moveUp2.setEnabled(false);
					moveDown2.setEnabled(false);
					top2.setEnabled(false);
					bottom2.setEnabled(false);
				}
			}
		};
		
		
		final JButton bottom=new ImageButton(SwingUtil.getImageGifIcon("moveBottom"));
		final JButton top = new ImageButton(SwingUtil.getImageGifIcon("moveTop"));
		
		final JButton moveUp = new JButton(SwingUtil.getImageGifIcon("upArrow"));
		moveUp.setToolTipText("<html>Click to move a used "+itemName+" up the " + leftListTip + " list</html>");
		final JButton moveDown=new JButton(SwingUtil.getImageGifIcon("downArrow"));
		moveDown.setToolTipText("<html>Click to move a used "+itemName+" down the " + leftListTip + " list</html>");
		final AbstractAction unselectAction = new AbstractAction() {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@SuppressWarnings("unchecked")
			public void actionPerformed(final ActionEvent e) {
				@SuppressWarnings("deprecation")
				final Object[] value = leftList.getSelectedValues();
				if (emptyNotAllowed && value.length == leftList.getModel().getSize()) {
					JOptionPane.showMessageDialog(dlg,"<html>Can not remove, as you need to have at least one " + itemName + " for column/tree display.</html>");
					return;
				}
				for(int i=0;i<value.length;i++){
					((DefaultListModel<T>)leftList.getModel()).removeElement(value[i]);
					((DefaultListModel<T>)rightList.getModel()).addElement((T)value[i]);
				}
				sort(((DefaultListModel<T>)rightList.getModel()));
				int[] indices=leftList.getSelectedIndices();

				for(int i=0;i<indices.length;i++){
					if(indices[i]==0){
						moveUp.setEnabled(false);
						top.setEnabled(false);
					}
					if(indices[i]==leftList.getModel().getSize()-1){
						moveDown.setEnabled(false);
						bottom.setEnabled(false);
					}
				}
				enableUpdates();
			}
		};

		rightButton = new JButton(unselectAction);
		rightButton.setToolTipText("<html>Click to <b><i>stop using</i></b> a selected "+itemNamePlural+" in "+leftListTip +" list </html>");
		final ListSelectionListener selectAndOrderListSelectionListener = new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent evt) {
				if (leftList.getSelectedIndex() > -1) {
					rightButton.setEnabled(true);
					ToolTipOnDemand.getSingleton().showLater(rightButton, false, null, -35, 45, false, 
							"<html><center><u>Double click</u> to <b><i>stop using</i></b><br><small>(or <u>drag</u> to "+rightListTip+" list<br>or <u>click</u> right button)</small></center><</html>");
				} else {
					rightButton.setEnabled(false);
				}
				top.setEnabled(true);
				bottom.setEnabled(true);
				moveUp.setEnabled(true);
				moveDown.setEnabled(true);
				int[] indices=leftList.getSelectedIndices();
				for(int i=0;i<indices.length;i++){
					if(indices[i]==0){
						top.setEnabled(false);
						moveUp.setEnabled(false);
					}
					else if(indices[i]==leftList.getModel().getSize()-1){
						bottom.setEnabled(false);
						moveDown.setEnabled(false);
					}
				}
				if(leftList.getModel().getSize() <=1 || indices.length <= 0) {
					top.setEnabled(false);
					bottom.setEnabled(false);
					moveUp.setEnabled(false);
					moveDown.setEnabled(false);
				}
				else if(leftList.getSelectedIndices().length == leftList.getModel().getSize()) {
					top.setEnabled(false);
					bottom.setEnabled(false);
					moveUp.setEnabled(false);
					moveDown.setEnabled(false);
				}
			}
		};
		
		rightButton.setVisible(this.isLeftRightOk);


		final AbstractAction selectAction = new AbstractAction() {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@SuppressWarnings("unchecked")
			public void actionPerformed(ActionEvent e) {
				@SuppressWarnings("deprecation")
				final Object[] value = rightList.getSelectedValues();
				for(int i=0;i<value.length;i++){
					((DefaultListModel<T>)leftList.getModel()).addElement((T)value[i]);
					((DefaultListModel<T>)rightList.getModel()).removeElement(value[i]);
					selectAndOrderListSelectionListener.valueChanged(new ListSelectionEvent(leftList,0,0,false));
				}
				int[] indices=leftList.getSelectedIndices();

				for(int i=0;i<indices.length;i++){
					if(indices[i]==0){
						top.setEnabled(false);
						moveUp.setEnabled(false);
					}
					if(indices[i]==leftList.getModel().getSize()-1){
						bottom.setEnabled(false);
						moveDown.setEnabled(false);
					}
				}
				if (ms!=null) {
					System.out.println("Marker sort now");
					sort(((DefaultListModel<T>)leftList.getModel()));
				}
				
				enableUpdates();
			}
		};


		anAction = new ActionListener() {
			public void actionPerformed(final ActionEvent a) {
				Mover.move(leftList, Mover.Direction.TOP);
				leftList.requestFocus();
				enableUpdates();
			}
		};
		SwingUtil.echoAction(leftList, top, anAction, 
				KeyStroke.getKeyStroke(KeyEvent.VK_HOME, InputEvent.ALT_MASK), 't');
		anAction = new ActionListener() {
			public void actionPerformed(final ActionEvent a) {
				Mover.move(leftList, Mover.Direction.BOTTOM);
				leftList.requestFocus();
				enableUpdates();
			}
		};
		SwingUtil.echoAction(leftList, bottom, anAction, 
				KeyStroke.getKeyStroke(KeyEvent.VK_END, InputEvent.ALT_MASK), 'b');
		moveUp.setMargin(new Insets(1, 2, 1, 2));
		moveDown.setMargin(new Insets(1, 2, 1, 2));
		anAction = new ActionListener() {
			public void actionPerformed(final ActionEvent a) {
				Mover.move(leftList, Mover.Direction.UP);
				leftList.requestFocus();
				enableUpdates();
			}
		};
		SwingUtil.echoAction(leftList, moveUp, anAction, 
				KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.ALT_MASK), 'u');
		anAction = new ActionListener() {
			public void actionPerformed(final ActionEvent a) {
				Mover.move(leftList, Mover.Direction.DOWN);
				leftList.requestFocus();
				enableUpdates();
			}
		};
		SwingUtil.echoAction(leftList, moveDown, anAction, 
				KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.ALT_MASK), 'd');
		final JPanel moveButtonPanel = new JPanel(new BorderLayout());
		final JPanel northPanel = new JPanel(new BorderLayout(0, 4));
		northPanel.add(top, "North");
		northPanel.add(moveUp, "South");
		moveButtonPanel.add(northPanel, "North");
		final JPanel southPanel = new JPanel(new BorderLayout(0, 4));
		southPanel.add(moveDown, "North");
		southPanel.add(bottom, "South");
		moveButtonPanel.add(southPanel, "South");
		mainPanel = new JPanel(new BorderLayout(0,0));
		mainPanel.setOpaque(true);
		rightListScrollPane = new JScrollPane(	rightList);
		rightListScrollPane.setPreferredSize(new Dimension(350, 250));

		final JPanel panel = new JPanel();
		panel.setLayout(new BorderLayout(0,0));
		panel.add(rightListScrollPane, BorderLayout.CENTER);
		panel.add(moveButtonPanel2, BorderLayout.EAST);
		
		final JPanel east=new JPanel(new BorderLayout(8, 1));
		east.setBorder(BorderFactory.createEmptyBorder(1, 2, 10, 2));
		final JPanel northEast=new JPanel();
		northEast.add(rightListHeader);
		east.add(northEast, BorderLayout.NORTH);
		east.add(panel, BorderLayout.CENTER);

		final JPanel actionPanel = new JPanel();
		actionPanel.setLayout(new BoxLayout(actionPanel, BoxLayout.Y_AXIS));

		leftButton = new JButton(selectAction);
		Basics.HearEnterKey(leftList, rightButton);
		Basics.HearEnterKey(rightList, leftButton);

		leftButton.setToolTipText("<html>Click to <b><i>use</i></b> selected "+itemNamePlural+" in list "+rightListTip);
		if (rightList != null) {
			rightList.addMouseListener(new MouseAdapter() {
				public void mouseClicked(final MouseEvent evt) {
					if (isLeftRightOk) {
						if (evt.getClickCount() == 2) { // Double-click
							selectAction.actionPerformed(new ActionEvent(rightList,
									evt.getID(), "select"));
						}
					}
				}
			});				
		}
		leftButton.setIcon(SwingUtil.getImageGifIcon("leftArrow"));
		if (rightList != null) {
			rightList.addListSelectionListener(selectAndOrderListSelectionListener2);
		}
		leftButton.setMargin(new Insets(1, 2, 1, 2));
		actionPanel.add(leftButton);
		leftList.setVisibleRowCount(rows);
		leftList.setDragEnabled(true);
		leftList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		leftButton.setEnabled(false);
		if (leftList != null) {
			if (doubleClickLeftList != null){
				leftList.setToolTipText(
						"<html><h3>More options exist ...</h3>"+
								"Double click on any item in the "+rightListTip+
						"<br>in order to access more options...</html>");
			} else {
				leftList.setToolTipText(
						"<html>To <b><i>stop using</i></b> a "+itemName+" <u>drag</u> it right to the "+rightListTip+
						" list or click right button.<br>To <b><i>reorder</i></b> a "+itemName+" drag up/down or click up/down buttons.</html>");
				leftLabel.setToolTipText(leftList.getToolTipText());
			}

			leftList.addMouseListener(new MouseAdapter() {
				public void mouseEntered(final MouseEvent evt){
					if (doubleClickLeftList != null){
						ToolTipManager.sharedInstance().setEnabled(false);
						ToolTipOnDemand.getSingleton().showWithoutCancelButton(
								leftList, false, evt.getX()+5, evt.getY()+5);							
					}
				}

				public void mouseExited(final MouseEvent evt){
					if (doubleClickLeftList != null){
						ToolTipManager.sharedInstance().setEnabled(true);
						ToolTipOnDemand.getSingleton().hideTipWindow();							
					}
				}

				public void mouseClicked(MouseEvent evt) {
					if (doubleClickLeftList != null){
						ToolTipManager.sharedInstance().setEnabled(true);
						ToolTipOnDemand.getSingleton().hideTipWindow();							
					}
					if (evt.getClickCount() == 2) { // Double-click
						if (isLeftRightOk) {
							if (doubleClickLeftList ==null){
								unselectAction.actionPerformed(new ActionEvent(leftList,
										evt.getID(), "select"));
							} else {
								final int index=leftList.locationToIndex(evt.getPoint());
								if (index>=0){
									doubleClickLeftList.respond(leftList.getModel().getElementAt(index), leftList);
									selectAndOrderListSelectionListener.valueChanged(new ListSelectionEvent(leftList,0,0,false));
								}
							}
						}
					}						
				}
			});				
		}		
		rightButton.setIcon(SwingUtil.getImageGifIcon("rightArrow"));
		if (leftList != null) {
			leftList
			.addListSelectionListener(selectAndOrderListSelectionListener);
		}
		rightButton.setEnabled(false);
		rightButton.setMargin(new Insets(1, 2, 1, 2));
		actionPanel.add(rightButton);
		leftListScrollPane= new JScrollPane(leftList);
		leftListScrollPane.setPreferredSize(new Dimension(350, 250));
		JPanel panel1 = new JPanel();
		panel1.setLayout(new BorderLayout(0,0));
		panel1.add(leftListScrollPane, BorderLayout.CENTER);
		panel1.add(moveButtonPanel, BorderLayout.EAST);
		final JPanel west=new JPanel(new BorderLayout(8,1));
		west.setBorder(BorderFactory.createEmptyBorder(1, 15, 10, 2));
		final JPanel northWest=new JPanel();
		if (leftLabel != null){
			northWest.add(leftLabel);
		} else {
			northWest.add(new JLabel("Selected"));
		}
		west.add(northWest, BorderLayout.NORTH);
		west.add(panel1, BorderLayout.CENTER);
		Box listsPanel = new Box(BoxLayout.LINE_AXIS);
		listsPanel.add(west);
		listsPanel.add(actionPanel);
		listsPanel.add(east);

		// Add the lists panel to the overall panel	
		if (guide!=null){
			final JPanel jp = new JPanel(new BorderLayout());
			jp.add(rightLabel, BorderLayout.NORTH);
			jp.add(new JPanel(), BorderLayout.CENTER);
			JLabel guideLbl = new JLabel(guide);
			jp.add(guideLbl, BorderLayout.SOUTH);
			mainPanel.add(jp, BorderLayout.NORTH);
		}
		mainPanel.add(listsPanel, BorderLayout.CENTER);
		final JPanel south=new JPanel(new BorderLayout());
		if (atBottom != null){
			buttonPanel.add(atBottom);
		}
		south.add(buttonPanel, BorderLayout.WEST);
		final JButton done=new JButton("Done");
		
		if (wnd != null){
			if (wnd instanceof JDialog) {
				final ActionListener taskPerformer=new ActionListener() {
					public void actionPerformed(ActionEvent evt) {
						SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								((JDialog)wnd).getRootPane().setDefaultButton(done);
								done.requestFocus();
								System.out.println("Focus is done...");
							}
						});
					}
				};
				javax.swing.Timer timer= new javax.swing.Timer(300 , taskPerformer);
				timer.setRepeats(false);
				timer.start();
			} 
		
			done.addActionListener(CpuInfo.getCloseAction(wnd));
			if (doneListener != null) {
				done.addActionListener(doneListener);
			}
			if (wnd instanceof RootPaneContainer){
				CpuInfo.registerEscape((RootPaneContainer)wnd, done);
			}
			if (updateButton!=null && showUpdate){
				//buttonPanel.add(updateButton); //Disabled as requested by Steve. W
				buttonPanel.add(cbAutoUpdate);
			}
			if (resetRequired) {
				final JButton resetButton=new JButton("Revert to original");
				resetButton.addActionListener(new ActionListener() {
					public void actionPerformed(final ActionEvent e) {
						leftModel.removeAllElements();
						rightModel.removeAllElements();
						Iterator<T> it=leftSet.iterator();
						while(it.hasNext()){
							leftModel.addElement(it.next());
						}
						it=rightSet.iterator();
						while(it.hasNext()){
							rightModel.addElement(it.next());
						}
						enableUpdates();
					}
				});
				buttonPanel.add(resetButton);
			}
			final JPanel doneJp=new JPanel();
			doneJp.add(done);
			south.add(doneJp, BorderLayout.EAST);
		}
		mainPanel.add(south, BorderLayout.SOUTH);
		if (leftList != null && leftList.getSelectedIndex() < 0) {
			top.setEnabled(false);
			bottom.setEnabled(false);
			moveUp.setEnabled(false);
			moveDown.setEnabled(false);
		}
		if (rightList != null && rightList.getSelectedIndex() < 0) {
			top2.setEnabled(false);
			bottom2.setEnabled(false);
			moveUp2.setEnabled(false);
			moveDown2.setEnabled(false);
		}
		setIsLeftRightOk(isLeftRightOk);
		btns=buttonPanel;
		return mainPanel;
	}
	public JPanel btns;
	ActionListener doneListener = null;

	public void setDoneListener(ActionListener doneListener) {
		this.doneListener = doneListener;
	}

	public void show(final boolean modal, final String location){
		if (location != null) {
			SwingUtil.position(owner, dlg, location);
		}
		dlg.setModal(modal);
		if (updateButton !=null){
			dlg.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					if (updateButton .isEnabled()){
						try{
							updateButton .doClick();
						}catch(final Exception ex){
							ex.printStackTrace(System.err);
						}
					}
					dlg.dispose();
				}
			});
		}
		dlg.setVisible(true);
		ToolTipManager.sharedInstance().setEnabled(true);		
	}

	public java.util.List<T> getResults(){
		final java.util.List<T> results = new ArrayList<T> ();
		for(int i=0;i<leftList.getModel().getSize();i++) {
			results.add((T)leftList.getModel().getElementAt(i));
		}
		return results;
	}

	boolean sortingAlphaCD=false;
	MarkerSorter ms;
	
	public void sortAlphaCD() {
		if (ms==null) {
			ms=new MarkerSorter();
		}
		resetLists((Collection) ms.sort(getLeftList()), (Collection) ms.sort(getRightList()));
	}
	
	public java.util.List<T> getLeftList(){
		final java.util.List<T> results = new ArrayList<T> ();
		for(int i=0;i<leftList.getModel().getSize();i++) {
			results.add((T)leftList.getModel().getElementAt(i));
		}
		return results;
	}

	public java.util.List<T> getRightList(){
		final java.util.List<T> results = new ArrayList<T> ();
		for(int i=0;i<rightList.getModel().getSize();i++) {
			results.add((T)rightList.getModel().getElementAt(i));
		}
		return results;
	}

	public String title="Select and order", guide="Pick the order of items you want";
	public void setLists(
			final Window owner,
			final Collection<T> leftItems,
			final Collection<T> rightItems,
			final String leftLabel,
			final String rightLabel) {
		final JLabel l1=SwingUtil.TitleLabel(new JLabel(leftLabel)), 
				l2=SwingUtil.TitleLabel(new JLabel(rightLabel));
		setLists(owner, leftItems, 
				rightItems, null,l1, l2, null, null, null);
	}
	
	public void setLists(
			final Window owner,
			final Collection<T> leftItems,
			final Collection<T> rightItems,
			final String leftLabel,
			final JComponent rightLabels) {
		final JLabel l1=SwingUtil.TitleLabel(new JLabel(leftLabel));
		setLists(owner, leftItems, 
				rightItems, null,l1, rightLabels, null, null, null);
	}
	
	public JPanel setLists(
			final Window owner,
			final Collection<T> leftItems,
			final Collection<T> rightItems,
			final String leftLabel,
			final JComponent rightLabels,
			final ListCellRenderer<T> leftLcr,
			final ListCellRenderer<T> rightLcr) {
		final JLabel l1=SwingUtil.TitleLabel(new JLabel(leftLabel));
		return setLists(owner, leftItems, 
				rightItems, null,l1, rightLabels, leftLcr, rightLcr, null);
	}
	
  	public static Collection<String> toCollection( final String arg1, final String ... args){
  		final Collection<String> c=new ArrayList<String>();
  		c.add(arg1);
  		for (final String arg:args){
  			c.add(arg);
  		}
  		return c;
  	}

	public static void main(String []args){
		final JButton u=new JButton("Show result");
		u.setToolTipText("Click to update results");
		
		final Orderer<String> or=new Orderer<String>(u, true);
		or.visibleRows=6;
		u.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				final java.util.List<String> c=or.getResults();
				System.out.println(c.size()+" selected..."+c);
			}
		});
		final Collection<String>rightItems=toCollection("File name (in file system)", "SPECIMEN", "$SRC", "TYPE", 
				"TUBE NAME", "keyword5", "keyword6", "k7", "k8", 
				"yada yaada yada");
		Collection<String> leftItems=toCollection("Sample keyword", "$SRC");
		or.title="Testing";
		or.setLists(
				null,  
				leftItems, 
				rightItems,
				"Used for name",
				"Available for name");
		or.btns.add(u);
		//or.setIsLeftRightOk(false);
		or.show(false, "south west");
	}
}
