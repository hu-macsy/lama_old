/*
 * MyGrappaListener.java
 *
 * Frame that displays the call graph.
 *
 * Created: 2006-02-20 Thomas Brandes <thomas.brandes@scai.fraunhofer.de>
 * Changed:
 *
 * $Id$
 *
 * Copyright (C) 2006 Fraunhofer SCAI, Germany
 *
 * All rights reserved
 *
 * http://www.scai.fhg.de/EP-CACHE/adaptor
 */

package adaptor.Calltree;

import java.util.List;
import java.util.Vector;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;

import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;

import org.apache.log4j.Logger;

import att.grappa.Element;
import att.grappa.GrappaBox;
import att.grappa.GrappaConstants;
import att.grappa.GrappaListener;
import att.grappa.GrappaPanel;
import att.grappa.GrappaPoint;
import att.grappa.Node;
import att.grappa.Subgraph;

/**
 * This class implements listeners on the grappa panel used for CalltreePM.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class MyGrappaListener implements GrappaConstants, GrappaListener, ActionListener
{

    /**
     * Common rror message used for problems on selections.
     */
    private static final String SELECTION_ERROR = "currentSelection improperly maintained";

    /**
     * This factor is used for zooming in and should be greater than 1.0.
     */
    private static final double ZOOM_IN_FACTOR = 1.25;

    /**
     * This factor is used for zooming out and should be less than 1.0.
     */
    private static final double ZOOM_OUT_FACTOR = 0.8;

    /**
     * These are the entries of the popup menu that appears with right
     * mouse selection within the grappa display.
     */
    private static final String[] MENU_ITEMS = {"Info", "Show", "Unfold",
                                  "Fold", "Call Environment", "Hide", "Zoom In", "Zoom Out", "Reset Zoom",
                                  "Scale to Fit"
                                               };

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( MyGrappaListener.class );

    /**
     * Pointer to the main class that implements the routines called for the actions.
     */
    private CTInterface myMain;

    /**
     * Constructor for our listener to the Grappa Panel. It needs the pointer
     * back to the main GUI that implements action routines.
     *
     * @param theMain is the pointer back to Calltree functions
     */
    public MyGrappaListener( CTInterface theMain )
    {

        this.myMain = theMain;

    }

    /**
     * This routine is used to hide a set of CT nodes. This
     * is done by marking them as false.
     *
     * @param nodes is an array with the nodes that will be hided.
     */
    private void hide( CTNode[] nodes )
    {

        if ( nodes != null )
        {

            // unmark all nodes

            for ( int i = 0; i < nodes.length; i++ )
            {

                nodes[i].setMarked( false );

            }

            if ( nodes.length > 0 )
            {

                myMain.redraw();
            }
        }

    }

    /**
     * This routine is used to unfold a set of CT nodes.
     *
     * @param nodes is an array of CT nodes that will be unfolded.
     */
    private void unfold( CTNode[] nodes )
    {

        if ( nodes != null )
        {

            // unfold all nodes of the array

            for ( int i = 0; i < nodes.length; i++ )
            {

                nodes[i].unfold();

            }

            if ( nodes.length > 0 )
            {

                myMain.redraw();
            }
        }

    }

    /**
     * Folding back a number of CT nodes.
     *
     * @param nodes us the array of calltree nodes that will be unfolded.
     */
    private void fold( CTNode[] nodes )
    {

        if ( nodes == null )
        {

            // no selection, so we fold all nodes

            logger.info( "reset all unfolded nodes" );

            myMain.resetUnfolding();

        }
        else
        {

            int count = 0;

            // fold back all selected nodes

            for ( int i = 0; i < nodes.length; i++ )
            {

                if ( nodes[i].fold() )
                {

                    count++;
                }

            }

            logger.info( count + " nodes were folded back" );

            if ( count > 0 )
            {

                myMain.redraw();
            }
        }

    }

    /**
     * This routine is called when an item of the popup menu
     * has been selected.
     *
     * @param item is the selected item
     * @param gp is the grappa panel where the selection was
     */
    public void doAction( JMenuItem item, GrappaPanel gp )
    {

        String cmd = ( ( JMenuItem ) item ).getText();

        if ( cmd.startsWith( "Reset" ) )
        {

            gp.setScaleToFit( false );
            gp.setScaleToSize( null );
            gp.resetZoom();
            gp.clearOutline();

        }
        else if ( cmd.startsWith( "Unfold" ) )
        {

            Subgraph subg = gp.getSubgraph();

            CTNode[] nodes = MySubgraph.getSelectedNodes( subg );

            unfold( nodes );


        }
        else if ( cmd.startsWith( "Hide" ) )
        {

            Subgraph subg = gp.getSubgraph();

            CTNode[] nodes = MySubgraph.getSelectedNodes( subg );

            hide( nodes );

        }
        else if ( cmd.startsWith( "Fold" ) )
        {

            Subgraph subg = gp.getSubgraph();

            CTNode[] nodes = MySubgraph.getSelectedNodes( subg );

            fold( nodes );

        }
        else if ( cmd.startsWith( "Call Environment" ) )
        {

            logger.info( "selected nodes + call environment" );

            myMain.showCallEnvironment();

        }
        else if ( cmd.startsWith( "Show" ) )
        {

            logger.info( "show the current selection" );

            Subgraph subg = gp.getSubgraph();

            CTNode node = MySubgraph.getSelectedNode( subg );

            if ( node != null )
            {

                node.showNode();
            }

        }
        else if ( cmd.startsWith( "Info" ) )
        {

            logger.info( "show the current selection" );

            Subgraph subg = gp.getSubgraph();

            CTNode node = MySubgraph.getSelectedNode( subg );

            myMain.showInfo( node );

        }
        else if ( cmd.startsWith( "Scale" ) )
        {

            gp.setScaleToFit( true );

        }
        else if ( cmd.startsWith( "Zoom In" ) )
        {

            gp.setScaleToFit( false );
            gp.setScaleToSize( null );
            gp.multiplyScaleFactor( ZOOM_IN_FACTOR );
            gp.clearOutline();

        }
        else if ( cmd.startsWith( "Zoom Out" ) )
        {

            gp.setScaleToFit( false );
            gp.setScaleToSize( null );
            gp.multiplyScaleFactor( ZOOM_OUT_FACTOR );
            gp.clearOutline();

        }

        logger.info( "doAction for " + cmd );

    } // doAction

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed( ActionEvent e )
    {

        Object actionSource = ( Object ) e.getSource();

        if ( actionSource instanceof JMenuItem )
        {

            Object parent = ( ( JMenuItem ) actionSource ).getParent();

            if ( parent instanceof JPopupMenu )
            {

                Object invoker = ( ( JPopupMenu ) parent ).getInvoker();

                if ( invoker instanceof GrappaPanel )
                {

                    // that is the constallation for an action

                    GrappaPanel gp = ( GrappaPanel ) invoker;

                    doAction( ( JMenuItem ) actionSource, gp );

                }

            }

        }

    } // actionPerformed

    /**
     * {@inheritDoc}
     *
     * @see att.grappa.GrappaListener#grappaTip(att.grappa.Subgraph, att.grappa.Element,
     *                                          att.grappa.GrappaPoint, int, att.grappa.GrappaPanel)
     */
    public String grappaTip( Subgraph subg, Element elem, GrappaPoint pt, int modifiers, GrappaPanel panel )
    {

        String tip = "";

        if ( elem == null )
        {

            tip = "Nothing selected";

        }
        else
        {

            switch ( elem.getType() )
            {

                case SUBGRAPH:
                    tip = "Subgraph";
                    break;
                case EDGE:
                    tip = "Edge: ";
                    break;
                case NODE:
                    Node node = ( Node ) elem;
                    CTNode ctnode = CalltreePM.getNode( node.getName() );
                    tip = "Node: " + node.getName() + " : " + ctnode.getCallNode().getNodeKindString();
                    break;
                default:
                    throw new RuntimeException( "unexpected type (" + elem.getType() + ")" );
            }
        }

        return tip;

    } // grappaTip

    /**
     * {@inheritDoc}
     */
    public void grappaDragged( Subgraph subg, GrappaPoint currentPt, int currentModifiers, Element pressedElem, GrappaPoint pressedPt,
                               int pressedModifiers, GrappaBox outline, GrappaPanel panel )
    {

        logger.debug( "MyGrappaListener: grappaDragged" );

    } // grappaDragged

    /**
     * {@inheritDoc}
     *
     * This routine is not used here.
     */
    public void grappaReleased( Subgraph subg, Element elem, GrappaPoint pt, int modifiers, Element pressedElem, GrappaPoint pressedPt,
                                int pressedModifiers, GrappaBox outline, GrappaPanel panel )
    {

        // logger.info("MyGrappaListener: grappaReleased");

    } // grapppaReleased

    /**
     * {@inheritDoc}
     *
     * @see att.grappa.GrappaListener#grappaPressed(att.grappa.Subgraph,
     * att.grappa.Element, att.grappa.GrappaPoint, int, att.grappa.GrappaPanel)
     */
    public void grappaPressed( Subgraph subg, Element elem, GrappaPoint pt, int modifiers, GrappaPanel panel )
    {

        int flag = modifiers & ( InputEvent.BUTTON2_MASK | InputEvent.BUTTON3_MASK );

        if ( flag != 0 && flag == modifiers )
        {

            // pop-up menu if button2 or button3

            JPopupMenu popup = new JPopupMenu();

            for ( int i = 0; i < MENU_ITEMS.length; i++ )
            {

                JMenuItem item = new JMenuItem( MENU_ITEMS[i] );
                item.addActionListener( this );
                popup.add( item );

            }

            java.awt.geom.Point2D mpt = panel.getTransform().transform( pt, null );
            popup.show( panel, ( int ) mpt.getX(), ( int ) mpt.getY() );

        } // mouse selection

    } // grappaPressed

    /**
     * This routine resets the selection in a subgraph.
     *
     * @param subg is the subgraph
     */
    private static void resetSelection( Subgraph subg )
    {

        if ( subg.currentSelection != null )
        {

            if ( subg.currentSelection instanceof Element )
            {

                ( ( Element ) ( subg.currentSelection ) ).highlight &= ~HIGHLIGHT_MASK;

            }
            else
            {

                List vec = ( Vector ) subg.currentSelection;

                for ( int i = 0; i < vec.size(); i++ )
                {

                    ( ( Element ) ( vec.get( i ) ) ).highlight &= ~HIGHLIGHT_MASK;

                } // for
            }

            subg.currentSelection = null;

        }
    }

    /**
     * This routine removes an element from the selected elements
     * of a subgraph.
     *
     * @param subg is the subgraph with a selection
     * @param selElem is the selected element
     */
    void removeSelectedElement( Subgraph subg, Element selElem )
    {

        // unselect selElem

        selElem.highlight &= ~SELECTION_MASK;

        if ( subg.currentSelection == null )
        {

            // something got messed up somewhere

            throw new InternalError( SELECTION_ERROR );

        }
        else if ( subg.currentSelection instanceof Element )
        {

            if ( ( ( Element ) ( subg.currentSelection ) ) != selElem )
            {

                // something got messed up somewhere

                throw new InternalError( SELECTION_ERROR );
            }

            subg.currentSelection = null;

        }
        else
        {

            // current selection has many elements, search for selElem

            List vec = ( Vector ) subg.currentSelection;

            boolean problem = true;

            for ( int i = 0; i < vec.size(); i++ )
            {

                if ( ( ( Element ) ( vec.get( i ) ) ) == selElem )
                {

                    vec.remove( i );
                    problem = false;
                    break;
                }

            } // for all elements of the vector

            if ( problem )
            {

                // something got messed up somewhere

                throw new InternalError( SELECTION_ERROR );
            }
        }
    }

    /**
     * This routine adds an element to the selection of a subgraph.
     *
     * @param subg is the subgraph
     * @param selElem is the selected element
     */
    private void addSelectedElement( Subgraph subg, Element selElem )
    {

        // select element

        selElem.highlight |= SELECTION_MASK;

        if ( subg.currentSelection == null )
        {

            subg.currentSelection = selElem;

        }
        else if ( subg.currentSelection instanceof Element )
        {

            Object obj = subg.currentSelection;
            subg.currentSelection = new Vector();
            ( ( Vector ) ( subg.currentSelection ) ).add( obj );
            ( ( Vector ) ( subg.currentSelection ) ).add( selElem );

        }
        else
        {

            ( ( Vector ) ( subg.currentSelection ) ).add( selElem );

        }
    }

    /**
     * {@inheritDoc}
     *
     * @see att.grappa.GrappaListener#grappaClicked(att.grappa.Subgraph,
     * att.grappa.Element, att.grappa.GrappaPoint, int, int, att.grappa.GrappaPanel)
     */
    public void grappaClicked( Subgraph subg, Element elem, GrappaPoint pt, int modifiers, int clickCount, GrappaPanel panel )
    {

        logger.info( "grappaClicked, clickcount = " + clickCount );

        boolean isSelection = ( modifiers == InputEvent.BUTTON1_MASK );

        boolean isExtension = ( modifiers == ( InputEvent.BUTTON1_MASK | InputEvent.CTRL_MASK ) );

        if ( !( isSelection || isExtension ) )
        {

            return;
        }

        Element selElem = elem;

        if ( selElem != null )
        {

            if ( selElem.getType() != NODE )
            {

                selElem = null;
            }
        }

        if ( selElem == null )
        {

            // so we have clicked the mouse somewhere

            if ( isSelection )
            {

                resetSelection( subg );
                subg.getGraph().repaint();

            }

        }
        else if ( selElem.getType() == NODE )
        {

            if ( clickCount == 2 )
            {

                // double click on a node, we open it

                Node node = ( Node ) selElem;

                logger.info( "Grappa Node " + node.getName() + " opened" );

            }
            else if ( isSelection )
            {

                // change selection only if it is a new element

                if ( !selElem.equals( subg.currentSelection ) )
                {

                    resetSelection( subg );

                    selElem.highlight |= SELECTION_MASK;
                    subg.currentSelection = selElem;
                    subg.getGraph().repaint();

                }

            }
            else if ( isExtension )
            {

                // Control Mouse-Click: remove or add element to the selection

                if ( ( selElem.highlight & SELECTION_MASK ) == SELECTION_MASK )
                {

                    removeSelectedElement( subg, selElem );

                }
                else
                {

                    addSelectedElement( subg, selElem );

                }

                subg.getGraph().repaint();
            }

        }

    } // grappaClicked

} // MyGrappaListener

