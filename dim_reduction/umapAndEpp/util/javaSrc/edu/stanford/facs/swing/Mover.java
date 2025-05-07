package edu.stanford.facs.swing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import javax.swing.JList;
import javax.swing.ListModel;

public class Mover {
	   public enum Direction {TOP, UP,DOWN, BOTTOM};
	   private static boolean contains(final int []a, final int arg) {
		   for (int i=0;i<a.length;i++) {
			   if (a[i]==arg) {
				   return true;
			   }
		   }
		   return false;
	   }
	   private static int []duplicate(final int[]in){
		   final int[]value=new int[in.length];
		   for (int i=0;i<value.length;i++) {
			   value[i]=in[i];
		   }
		   return value;
	   }
	   
	   public static int []move(final int []selections, final Direction direction, final int lastIndex){
		   final int []oldLocations=duplicate(selections);
		   Arrays.sort(oldLocations);
		   return _move(oldLocations, direction, lastIndex);
	   }
	   
	   private static int []_move(final int []oldLocations, final Direction direction, final int lastIndex){
		   final int []newLocations;
		   switch(direction) {
		   case TOP:
			   newLocations=moveTop(oldLocations, 0);
			   break;
		   case BOTTOM:
			   newLocations=moveBottom(oldLocations, lastIndex);
			   break;
		   case UP:
			   newLocations=moveUp(oldLocations, 0);
			   break;
		   default:
			   newLocations=moveDown(oldLocations, lastIndex);
			   break;
		   }
		   return newLocations;
	   }
	   
	   public static Object[]move(final Object []in, final int []selections, final Direction direction){
		   final int []oldLocations=duplicate(selections);
		   Arrays.sort(oldLocations);
		   final int []newLocations=_move(oldLocations, direction, in.length-1);
		   final Object[]value=move(in, oldLocations, newLocations);
		   return value;
	   }
	   
	   private static int []moveUp(final int []input, final int firstIndex) {
		   final int []s=duplicate(input);
		   // s is sorted
		   int f=firstIndex;
		   for (int i=0;i<s.length;i++) {
			   if (s[i] != f) {
				   break;
			   }
			   f++;
		   }
		   for (int i=0;i<s.length;i++) {
			   if (s[i] > f) {
				   s[i]--;
			   }
		   }
		   return s;
	   }

	   private static int []moveDown(final int []input, final int lastIndex) {
		   final int []s=duplicate(input);
		   int l=lastIndex;
		   for (int i=s.length-1;i>=0;i--) {
			   if (s[i] != l) {
				   break;
			   }
			   l--;
		   }
		   for (int i=s.length-1;i>=0;i--) {
			   if (s[i] < l) {
				   s[i]++;
			   }
		   }
		   return s;
	   }
	   private static int []moveBottom(final int []input, final int bottomIndex) {
		   final int []s=duplicate(input);
		   
		   int j=0;
		   for (int i=s.length-1;i>=0;i--) {
			    s[i]=bottomIndex-j;
			   j++;
		   }
		   return s;
	   }
	   private static int []moveTop(final int []input, final int topIndex) {
		   final int []s=duplicate(input);
		   for (int i=0;i<input.length;i++) {
			    s[i]=topIndex+i;
			   
		   }
		   return s;
	   }
	   
	   static Object []move(final Object[] o, final int []oldLocations, final int[]newLocations) {
		   final Object []v=new Object[o.length];
		   for (int i=0;i<oldLocations.length;i++) {
			   v[ newLocations[i]]=o[ oldLocations[i]];
		   }
		   int j=0;
		   for (int i=0;i<o.length;i++) {
			   if (!contains(oldLocations, i)) {
			   while (v[j] != null) {
				   j++;
			   }
			   v[j]=o[i];
			   }
		   }
		   return v;
	   }
	   
	    public static Object []getDataList(final JList jl) {
	 	   final ListModel lm=jl.getModel();
	 	   final Collection c=new ArrayList();
	 	   for (int i=0;i<lm.getSize();i++) {
	 		   c.add(lm.getElementAt(i));
	 	   }   
	 	   return c.toArray();
	    }

	    public static void move(final JList jl, Direction d) {
	    	final Object []v1=getDataList(jl);
	    	final int []selected=jl.getSelectedIndices();
	    	final Object []v2=move(v1, selected, d);
	    	final int []n=move(selected, d, v1.length-1);
	    	jl.clearSelection();
	    	jl.setListData(v2);
	    	for (final int i:n) {
	    		jl.addSelectionInterval(i, i);    		
	    	}
	    	if (d == Direction.UP || d== Direction.TOP) {
	    		jl.ensureIndexIsVisible(n[0]);
	    	}else {
	    		jl.ensureIndexIsVisible(n[n.length-1]);        		
	    	}
	    }

}
