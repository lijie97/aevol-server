//*****************************************************************************
// S.P.E.A.R. library - Simulator of Physical Environment for Animat Research
// Copyright (C) 2003  S�bastien GRIPON, Fran�ois PERRIN, H�di SOULA, Virginie
//                     MATHIVET - PRISMa
// Web: http://prisma.insa-lyon.fr
// Original Authors : Guillaume BESLON, H�di SOULA
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************

// File Information
// $Id: ae_list.h,v 1.6 2005/10/25 17:18:45 cknibbe Exp $

#ifndef __AE_LIST_H__
#define __AE_LIST_H__


#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <ae_object.h>



class ae_list;

enum ae_delete
{
  NO_DELETE  = 0,
  DELETE_OBJ = 1
};





class ae_list_node
{
  friend class ae_list;
  public :
    inline ae_list_node( ae_object* obj, ae_list_node* prev = NULL, ae_list_node* next = NULL );
  
    // Make a new node sharing the same object as 'model'
    // The created node is unlinked
    inline ae_list_node( const ae_list_node &model );
    
    virtual inline ~ae_list_node();

    inline ae_object*    get_obj( void );
    inline ae_list_node* get_prev( void );
    inline ae_list_node* get_next( void );
    
    inline void set_obj( ae_object* obj );

  protected :
    ae_list_node( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    
    ae_object*     _obj;
    ae_list_node*  _prev;
    ae_list_node*  _next;
};




class ae_list : public ae_object
{
  public :
    inline ae_list();
    //~ inline ae_list( ae_list* parent );
    inline ae_list( const ae_list &model );
    inline virtual ~ae_list();

    // Add obj in a newly created node at the end (resp beginning) of the list
    // => syntaxic sugar for add_after( obj, _last );
    // (resp add_before( obj, _first );)
    inline ae_list_node* add( ae_object* obj );
    inline ae_list_node* add_front( ae_object* obj );

    // Link a node at the end (resp beginning) of the list. The node must be unlinked
    // => syntaxic sugar for insert_after( node, _last );
    // (resp insert_before( obj, _first );)
    inline void add( ae_list_node* node );
    inline void add_front( ae_list_node* node );

    // Add obj in a newly created node before (resp after) node n
    inline ae_list_node* add_before( ae_object* obj, ae_list_node* n );
    inline ae_list_node* add_after(  ae_object* obj, ae_list_node* n );

    // Insert node n1 before (resp after) node n2
    inline void insert_before( ae_list_node * n1, ae_list_node * n2 ); // TODO add au lieu de insert
    inline void insert_after(  ae_list_node * n1, ae_list_node * n2 );

    // Invert nodes n1 and n2
    inline void invert( ae_list_node * n1, ae_list_node * n2 );

    // Remove object from list and delete it if delete_obj == DELETE_OBJ
    inline void remove( ae_object* obj,     ae_delete delete_node, ae_delete delete_obj );
    inline void remove( ae_list_node* node, ae_delete delete_node, ae_delete delete_obj );
    
    // Remove all the nodes between node_1 and node_2 (included) from the list and return
    // a new list containing the removed nodes (order is conserved)
    inline ae_list* extract_sublist( ae_list_node* first_node, ae_list_node * last_node );
    inline ae_list* extract_starting_sublist( ae_list_node* last_node );
    inline ae_list* extract_ending_sublist( ae_list_node* first_node );

    // Erase the list (nodes are deleted). Objects are deleted if delete_obj == DELETE_OBJ
    inline void erase( ae_delete delete_obj );

    // Add the elements of 'append' at the end of the list.
    // 'append' is emptied.
    inline void merge( ae_list* append );
    
    // Add the elements of 'to_add' as new nodes at the end of the list.
    // 'to_add' is untouched.
    inline void add_list( ae_list* const to_add );


    inline int32_t  get_nb_elts( void ) const;
    inline bool     is_empty( void ) const;

    inline int32_t    get_position( ae_object* obj ) const;
    inline ae_list_node*     get_node( int32_t pos ) const;
    inline ae_object* get_object( int32_t pos ) const;

    // Accessors
    inline ae_list_node* get_first( void ) { return _first; };
    inline ae_list_node* get_last( void ) { return _last; };

    // Search for a value WITHIN the object.
    inline ae_list_node* bsearch( void* needle, int ( * comparator ) ( const void * value, const void * object ) ) const;


  protected :
    ae_list_node* _first;
    ae_list_node* _last;
    int32_t       _nb_elts;
};


// Constructor
ae_list::ae_list()
{
  _last = _first = NULL;
  _nb_elts = 0;
}


// Doesn't create new objects
ae_list::ae_list( const ae_list &model )
{
  _last = _first = NULL;
  _nb_elts = 0;

  ae_list_node* node = model._first;

  while ( node != NULL )
  {
    add( node->get_obj() );

    node = node->get_next();
  }
}


// Destructor
ae_list::~ae_list()
{
  erase( NO_DELETE );
}

// Add obj in a newly created node at the end of the list
ae_list_node* ae_list::add( ae_object* obj )
{
  ae_list_node* node = new ae_list_node( obj, _last, NULL );
  assert( node );
  
  if ( _last != NULL ) _last->_next = node;
  _last = node;
  if ( _first == NULL ) _first = node;
  _nb_elts++;
  return node;
}

// Add obj in a newly created node at the beginning of the list
ae_list_node* ae_list::add_front( ae_object* obj )
{
  ae_list_node* node = new ae_list_node( obj, NULL, _first );
  assert( node );
  
  if ( _first != NULL ) _first->_prev = node;
  _first = node;
  if ( _last == NULL ) _last = node;
  _nb_elts++;
  return node;
}


// Add a node at the end of the list. The node must be unlinked
void ae_list::add( ae_list_node* node )
{
  node->_prev = _last;
  node->_next = NULL;
  if( _last != NULL ) _last->_next = node;
  _last = node;
  if( _first == NULL ) _first = node;
  _nb_elts++;
}


// Add a node at the beginning of the list. The node must be unlinked
void ae_list::add_front( ae_list_node* node )
{
  node->_prev = NULL;
  node->_next = _first;
  if( _first != NULL ) _first->_prev = node;
  _first = node;
  if( _last == NULL ) _last = node;
  _nb_elts++;
}


// Add obj in a newly created node before node n
ae_list_node* ae_list::add_before( ae_object* obj, ae_list_node* n )
{
  ae_list_node* node = new ae_list_node( obj, NULL, NULL );
  assert( node );
  
  insert_before( node, n );
  
  return node;
}


// Add obj in a newly created node after node n
ae_list_node* ae_list::add_after( ae_object* obj, ae_list_node* n )
{
  ae_list_node* node = new ae_list_node( obj, NULL, NULL );
  assert( node );
  
  insert_after( node, n );
  
  return node;
}


// Insert node n1 before n2
void ae_list::insert_before( ae_list_node * n1, ae_list_node * n2 )
{
  // save node before n2
  ae_list_node * n2_prev = n2->_prev;

  // link n1 -> n2
  n1->_next = n2;
  n2->_prev = n1;

  // link (n2->_prev) -> n1
  n1->_prev = n2_prev;
  if (n2_prev == NULL)
  {
    // n2 was _first
    _first = n1;
  }
  else
  {
    n2_prev->_next = n1;
  }

  _nb_elts++;
}


// Insert node n1 after n2
void ae_list::insert_after( ae_list_node * n1, ae_list_node * n2 )
{
  assert( n1 != NULL && n2 != NULL );
  
  // save node after n2
  ae_list_node * n2_next = n2->_next;

  // link n2 -> n1
  n2->_next = n1;
  n1->_prev = n2;

  // link n1 -> (n2->_next)
  n1->_next = n2_next;
  if (n2_next == NULL)
  {
    // n2 was _last
    _last = n1;
  }
  else
  {
    n2_next->_prev = n1;
  }

  _nb_elts++;
}


// Invert nodes n1 and n2
void ae_list::invert( ae_list_node* n1, ae_list_node* n2 )
{
  // Optim if nodes are consecutive
  if ( n2 == n1->_next || n1 == n2->_next )
  {
    if ( n1 == n2->_next )
    {
      // Exchange n1 and n2
      ae_list_node* tmp = n1;
      n1 = n2;
      n2 = tmp;
    }

    if ( n1 != _first ) n1->_prev->_next = n2;
    else _first = n2;
    if ( n2 != _last ) n2->_next->_prev = n1;
    else _last = n1;
    n2->_prev = n1->_prev;
    n1->_next = n2->_next;
    n2->_next = n1;
    n1->_prev = n2;
  }
  else
  {
     // TODO
    printf( "NOT IMPLEMENTED!!! %s %d\n", __FILE__, __LINE__ );
    getchar();
  }
}


// Remove object from list and delete it if delete_obj == DELETE_OBJ
void ae_list::remove( ae_object* obj, ae_delete delete_node, ae_delete delete_obj )
{
  //printf("trying to rmove an object, nb_el = %ld\n", _nb_elts);
  for ( ae_list_node* node = _first ; node != NULL ; node = node->_next )
  {
    //printf("scanning obj %p for obj %p\n",node->_obj, obj);
    if ( node->_obj == obj )
    {
      // unlink node
      if ( node == _first )
      {
        _first = node->_next;
      }
      else if ( node->_prev )
      {
        node->_prev->_next = node->_next;
      }

      if ( node == _last )
      {
        _last = node->_prev;
      }
      else if ( node->_next )
      {
        node->_next->_prev = node->_prev;
      }

      // delete object (if requested) and node
      if( delete_obj )
      {
        if ( node->_obj != NULL )
        {
          delete node->_obj;
          node->_obj = NULL;
        }
        else
        {
          printf( "possible attempt to delete something twice: %s %d\n", __FILE__, __LINE__ );
        }
      }

      if ( delete_node )
      {
        if ( node != NULL )
        {
          delete node;
          node = NULL;
        }
        else
        {
          printf( "possible attempt to delete something twice: %s %d\n", __FILE__, __LINE__ );
        }
      }

      _nb_elts--;
      return;
    }

  }

  printf("WARNING: object to remove not found in the list ! \n");
}


// Remove node from list and delete the corrresponding object if delete_obj == DELETE_OBJ
void ae_list::remove( ae_list_node* node, ae_delete delete_node, ae_delete delete_obj )
{
  // unlink node
  if( node == _first )
  {
    _first = node->_next;
  }
  else if(node->_prev) node->_prev->_next = node->_next;

  if(node == _last)
  {
    _last = node->_prev;
  }
  else if(node->_next) node->_next->_prev = node->_prev;

  // delete object (if requested) and node
  if ( delete_obj )
  {
    if ( node->_obj != NULL )
    {
      delete node->_obj;
      node->_obj = NULL;
    }
    else printf("possible attempt to delete something twice: %s %d\n", __FILE__, __LINE__);
  }

  if ( delete_node )
  {
    if ( node != NULL )
    {
      delete node;
      node = NULL;
    }
    else
    {
      printf("possible attempt to delete something twice: %s %d\n", __FILE__, __LINE__);
    }
  }

  _nb_elts--;
}

ae_list* ae_list::extract_sublist( ae_list_node* first_node, ae_list_node* last_node )
{
  assert( first_node != NULL && last_node != NULL );
  
  ae_list* new_list = new ae_list();
  
  new_list->_first  = first_node;
  new_list->_last   = last_node;
  
  if ( first_node == _first && last_node == _last )
  {
    new_list->_nb_elts  = _nb_elts;
    
    _first    = NULL;
    _last     = NULL;
    _nb_elts  = 0;
  }
  else
  {
    if ( first_node == _first )
    {
      _first = last_node->_next;
      _first->_prev = NULL;
      last_node->_next = NULL;
    }
    else if ( last_node == _last )
    {
      _last = first_node->_prev;
      _last->_next = NULL;
      first_node->_prev = NULL;
    }
    else
    {
      first_node->_prev->_next = last_node->_next;
      last_node->_next->_prev  = first_node->_prev;
      first_node->_prev = NULL;
      last_node->_next = NULL;
    }
    
    
    // Update number of elements of both lists
    ae_list_node* node = first_node;
    new_list->_nb_elts++;
    _nb_elts--;
    
    while ( node != last_node )
    {
      new_list->_nb_elts++;
      _nb_elts--;
      
      node = node->_next;
    } 
  }
  
  return new_list;
}

ae_list* ae_list::extract_starting_sublist( ae_list_node* last_node )
{
  ae_list* new_list = new ae_list();
  
  new_list->_first  = _first;
  new_list->_last   = last_node;
  
  if ( last_node == _last )
  {
    new_list->_nb_elts  = _nb_elts;
    
    _first    = NULL;
    _last     = NULL;
    _nb_elts  = 0;
  }
  else
  {
    _first = last_node->_next;
    
    // Update number of elements of both lists
    ae_list_node* node = _first;
    new_list->_nb_elts++;
    _nb_elts--;
    
    while ( node != last_node )
    {
      new_list->_nb_elts++;
      _nb_elts--;
      
      node = node->_next;
    } 
  }
  
  return new_list;
}

ae_list* ae_list::extract_ending_sublist( ae_list_node* first_node )
{
  ae_list* new_list = new ae_list();
  
  new_list->_first  = first_node;
  new_list->_last   = _last;
  
  if ( first_node == _first )
  {
    new_list->_nb_elts  = _nb_elts;
    
    _first    = NULL;
    _last     = NULL;
    _nb_elts  = 0;
  }
  else
  {
    _last = first_node->_prev;
    
    // Update number of elements of both lists
    ae_list_node* node = first_node;
    new_list->_nb_elts++;
    _nb_elts--;
    
    while ( node != NULL )
    {
      new_list->_nb_elts++;
      _nb_elts--;
      
      node = node->_next;
    } 
  }
  
  return new_list;
}

// Erase the list (nodes are deleted). Objects are deleted if delete_obj == DELETE_OBJ
void ae_list::erase( ae_delete delete_obj )
{
  //  printf("nb elts in list %d\n", _nb_elts);

  ae_list_node* node = _first;

  while ( node != NULL )
  {
    ae_list_node* next = node->_next;

    if( delete_obj )
    {
      if ( node->_obj != NULL )
      {
        delete node->_obj;
        node->_obj = NULL;
      }
      else
      {
        printf( "possible attempt to delete something twice: %s %d\n", __FILE__, __LINE__ );
      }
    }
    
    delete node;
    node = next;
  }

  _last = _first = NULL;
  _nb_elts = 0;
}

void ae_list::merge( ae_list* append )
{
  if ( append == NULL || append->is_empty() ) return; // Nothing to do

  if ( this->is_empty() )
  {
    // Copy 'append'
    _first    = append->_first;
    _last     = append->_last;
    _nb_elts  = append->_nb_elts;

    // Empty 'append'
    append->_first    = NULL;
    append->_last     = NULL;
    append->_nb_elts  = 0;
  }
  else
  {
    // Link the lists together
    _last->_next = append->_first;
    append->_first->_prev = _last;

    // Update '_last'
    _last = append->_last;

    // Update number of elements
    _nb_elts += append->_nb_elts;

    // Empty 'append'
    append->_first    = NULL;
    append->_last     = NULL;
    append->_nb_elts  = 0;
  }
}

inline void ae_list::add_list( ae_list* const to_add )
{
  if ( to_add == NULL ) return; // Nothing to do
  
  ae_list_node* to_add_node = to_add->get_first();
  
  while ( to_add_node != NULL )
  {
    add( new ae_list_node( *to_add_node ) );
    
    to_add_node = to_add_node->get_next();
  }
}

bool ae_list::is_empty( void ) const
{
  return ( _first == NULL );
}


int32_t ae_list::get_nb_elts( void ) const
{
  return _nb_elts;
}


/*!
  Returns the position of the object provided or -1 if not found
 */
int32_t ae_list::get_position( ae_object* obj ) const
{
  int32_t i = 0;

  ae_list_node * node = _first;
  while ( node != NULL )
  {
    if ( node->_obj == obj ) return i;
    i++;
    node = node->_next;
  }

  return -1;
}


/*!
  Returns the node at position <pos> or NULL if pos is invalid
 */
ae_list_node* ae_list::get_node( int32_t pos ) const
{
  if ( pos >= _nb_elts ) return NULL;

  ae_list_node * node = _first;
  for ( int32_t i = 0 ; i < pos ; i++ )
  {
    node = node->_next;
  }
  return node;
}


/*!
  Returns the object at position <pos> or NULL if pos is invalid
 */
ae_object* ae_list::get_object( int32_t pos ) const
{
  if ( pos >= _nb_elts ) return NULL;

  ae_list_node * node = _first;
  for ( int32_t i = 0 ; i < pos ; i++ )
  {
    node = node->_next;
  }
  return node->_obj;
}


ae_list_node* ae_list::bsearch( void* needle, int ( * comparator ) ( const void * value, const void * object ) ) const
// Returns a pointer to the first object in the list where needle was found
// Returns NULL if needle was not found
// comparator : function that returns true if <value> was found in <object>
{
  ae_list_node* node = _first;

  while ( node != NULL )
  {
    if ( (*comparator)( needle, (void*)node->get_obj() ) == 0 ) return node;

    node = node->_next;
  }

  return NULL;
}




//******************************************************************************
//                                  ae_list_node
//******************************************************************************
inline ae_list_node::ae_list_node( ae_object* obj, ae_list_node* prev, ae_list_node* next )
{
  _obj  = obj;
  _prev = prev;
  _next = next;
};

ae_list_node::ae_list_node( const ae_list_node &model )
{
  _obj  = model._obj;
  _prev = NULL;
  _next = NULL;
};

ae_list_node::~ae_list_node()
{
}

inline ae_object* ae_list_node::get_obj( void )
{
  return _obj;
}

inline ae_list_node* ae_list_node::get_prev( void )
{
  return _prev;
}

inline ae_list_node* ae_list_node::get_next( void )
{
  return _next;
}

inline void ae_list_node::set_obj( ae_object* obj )
{
  _obj = obj;
}

#endif // __AE_LIST_H__
