// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
// 
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 
// ****************************************************************************


#ifndef AE_LIST_H
#define AE_LIST_H


#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>



template <typename T>
class ae_list;




template <typename T>
class ae_list_node
{
  friend class ae_list<T>;
  
  public :
    inline ae_list_node( T obj, ae_list_node* prev = NULL, ae_list_node* next = NULL );
  
    // Make a new node sharing the same object as 'model'
    // The created node is unlinked
    inline ae_list_node( const ae_list_node &model );
    
    virtual inline ~ae_list_node( void );

    inline T&             get_obj( void );
    inline ae_list_node*  get_prev( void );
    inline ae_list_node*  get_next( void );
    
    inline void set_obj( T obj );

  protected :
    ae_list_node( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    
    T _obj;
    ae_list_node* _prev;
    ae_list_node* _next;
};




template <typename T>
class ae_list
{
  public :
    inline ae_list();
    //~ inline ae_list( ae_list* parent );
    /*inline ae_list( const ae_list &model );*/
    inline virtual ~ae_list();

    // Add obj in a newly created node at the end (resp beginning) of the list
    // => syntaxic sugar for add_after( obj, _last );
    // (resp add_before( obj, _first );)
    inline ae_list_node<T>* add( T obj );
    inline ae_list_node<T>* add_front( T obj );

    // Link a node at the end (resp beginning) of the list. The node must be unlinked
    // => syntaxic sugar for insert_after( node, _last );
    // (resp insert_before( obj, _first );)
    inline void add( ae_list_node<T>* node );
    inline void add_front( ae_list_node<T>* node );

    // Add obj in a newly created node before (resp after) node n
    inline ae_list_node<T>* add_before( T obj, ae_list_node<T>* n );
    inline ae_list_node<T>* add_after(  T obj, ae_list_node<T>* n );

    // Insert node n1 before (resp after) node n2
    inline void insert_before( ae_list_node<T>* n1, ae_list_node<T>* n2 ); // TODO add au lieu de insert
    inline void insert_after(  ae_list_node<T>* n1, ae_list_node<T>* n2 );

    // Invert nodes n1 and n2
    inline void invert( ae_list_node<T>* n1, ae_list_node<T>* n2 );

    // Remove object from list and delete it if delete_obj is true
    inline void remove( T* obj, bool delete_node, bool delete_obj );
    inline void remove( ae_list_node<T>* node, bool delete_node, bool delete_obj );
    
    // Remove all the nodes between node_1 and node_2 (included) from the list and return
    // a new list containing the removed nodes (order is conserved)
    inline ae_list<T> * extract_sublist( ae_list_node<T> * first_node, ae_list_node<T> * last_node );
    inline ae_list<T> * extract_starting_sublist( ae_list_node<T> * last_node );
    inline ae_list<T> * extract_ending_sublist( ae_list_node<T> * first_node );

    // Erase the list (nodes are deleted). Objects are deleted if delete_obj == DELETE_OBJ
    inline void erase( bool delete_obj );

    // Add the elements of 'append' at the end of the list.
    // 'append' is emptied.
    inline void merge( ae_list<T> * append );
    
    // Add the elements of 'to_add' as new nodes at the end of the list.
    // 'to_add' is untouched.
    inline void add_list( ae_list<T> * const to_add );


    inline int32_t  get_nb_elts( void ) const;
    inline bool     is_empty( void ) const;

    inline int32_t        get_position( T* obj ) const;
    inline ae_list_node<T> * get_node( int32_t pos ) const;
    inline T*             get_object( int32_t pos ) const;

    // Accessors
    inline ae_list_node<T>* get_first( void ) { return _first; };
    inline ae_list_node<T>* get_last( void )  { return _last; };

    // Search for a value WITHIN the object.
    inline ae_list_node<T>* bsearch( void* needle, int ( * comparator ) ( const void * value, const void * object ) ) const;


  protected :
    ae_list_node<T>* _first;
    ae_list_node<T>* _last;
    int32_t       _nb_elts;
};


// Constructor
template <typename T>
ae_list<T>::ae_list()
{
  _last = _first = NULL;
  _nb_elts = 0;
}


// Doesn't create new objects
/*ae_list::ae_list( const ae_list &model )
{
  _last = _first = NULL;
  _nb_elts = 0;

  ae_list_node* node = model._first;

  while ( node != NULL )
  {
    add( node->get_obj() );

    node = node->get_next();
  }
}*/


// Destructor
template <typename T>
ae_list<T>::~ae_list()
{
  erase( false );
}

// Add obj in a newly created node at the end of the list
template <typename T>
ae_list_node<T>* ae_list<T>::add( T obj )
{
  ae_list_node<T>* node = new ae_list_node<T>( obj, _last, NULL );
  assert( node );
  
  if ( _last != NULL ) _last->_next = node;
  _last = node;
  if ( _first == NULL ) _first = node;
  _nb_elts++;
  return node;
}

// Add obj in a newly created node at the beginning of the list
template <typename T>
ae_list_node<T>* ae_list<T>::add_front( T obj )
{
  ae_list_node<T>* node = new ae_list_node<T>( obj, NULL, _first );
  assert( node );
  
  if ( _first != NULL ) _first->_prev = node;
  _first = node;
  if ( _last == NULL ) _last = node;
  _nb_elts++;
  return node;
}


// Add a node at the end of the list. The node must be unlinked
template <typename T>
void ae_list<T>::add( ae_list_node<T>* node )
{
  node->_prev = _last;
  node->_next = NULL;
  if( _last != NULL ) _last->_next = node;
  _last = node;
  if( _first == NULL ) _first = node;
  _nb_elts++;
}


// Add a node at the beginning of the list. The node must be unlinked
template <typename T>
void ae_list<T>::add_front( ae_list_node<T>* node )
{
  node->_prev = NULL;
  node->_next = _first;
  if( _first != NULL ) _first->_prev = node;
  _first = node;
  if( _last == NULL ) _last = node;
  _nb_elts++;
}


// Add obj in a newly created node before node n
template <typename T>
ae_list_node<T>* ae_list<T>::add_before( T obj, ae_list_node<T>* n )
{
  ae_list_node<T>* node = new ae_list_node<T>( obj, NULL, NULL );
  assert( node );
  
  insert_before( node, n );
  
  return node;
}


// Add obj in a newly created node after node n
template <typename T>
ae_list_node<T>* ae_list<T>::add_after( T obj, ae_list_node<T>* n )
{
  ae_list_node<T>* node = new ae_list_node<T>( obj, NULL, NULL );
  assert( node );
  
  insert_after( node, n );
  
  return node;
}


// Insert node n1 before n2
template <typename T>
void ae_list<T>::insert_before( ae_list_node<T> * n1, ae_list_node<T> * n2 )
{
  // save node before n2
  ae_list_node<T> * n2_prev = n2->_prev;

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
template <typename T>
void ae_list<T>::insert_after( ae_list_node<T> * n1, ae_list_node<T> * n2 )
{
  assert( n1 != NULL && n2 != NULL );
  
  // save node after n2
  ae_list_node<T> * n2_next = n2->_next;

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
template <typename T>
void ae_list<T>::invert( ae_list_node<T> * n1, ae_list_node<T> * n2 )
{
  // Optim if nodes are consecutive
  if ( n2 == n1->_next || n1 == n2->_next )
  {
    if ( n1 == n2->_next )
    {
      // Exchange n1 and n2
      ae_list_node<T> * tmp = n1;
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


// Remove object from list and delete it if delete_obj is true
template <typename T>
void ae_list<T>::remove( T* obj, bool delete_node, bool delete_obj )
{
  //printf("trying to remove an object, nb_el = %ld\n", _nb_elts);
  for ( ae_list_node<T> * node = _first ; node != NULL ; node = node->_next )
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

      // delete object (if requested)
      if ( delete_obj )
      {
        if ( node->_obj != NULL )
        {
          delete node->_obj;
          node->_obj = NULL;
        }
        else
        {
          printf( "%s:%d: Warning: Possible attempt to delete something twice\n", __FILE__, __LINE__ );
        }
      }

      // delete node (if requested)
      if ( delete_node )
      {
        if ( node != NULL )
        {
          delete node;
          node = NULL;
        }
        else
        {
          printf( "%s:%d: Warning: Possible attempt to delete something twice\n", __FILE__, __LINE__ );
        }
      }

      _nb_elts--;
      return;
    }
  }
}


// Remove node from list and delete the corrresponding object if delete_obj is true
template <typename T>
void ae_list<T>::remove( ae_list_node<T> * node, bool delete_node, bool delete_obj )
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

  // delete object (if requested)
  if ( delete_obj )
  {
    if ( node->_obj != NULL )
    {
      delete node->_obj;
      node->_obj = NULL;
    }
    else
    {
      printf( "%s:%d: Warning: Possible attempt to delete something twice\n", __FILE__, __LINE__ );
    }
  }

  // delete node (if requested)
  if ( delete_node )
  {
    if ( node != NULL )
    {
      delete node;
      node = NULL;
    }
    else
    {
      printf( "%s:%d: Warning: Possible attempt to delete something twice\n", __FILE__, __LINE__ );
    }
  }

  _nb_elts--;
}

template <typename T>
ae_list<T> * ae_list<T>::extract_sublist( ae_list_node<T> * first_node, ae_list_node<T> * last_node )
{
  assert( first_node != NULL && last_node != NULL );
  
  ae_list<T> * new_list = new ae_list<T>();
  
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
    ae_list_node<T> * node = first_node;
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

template <typename T>
ae_list<T> * ae_list<T>::extract_starting_sublist( ae_list_node<T> * last_node )
{
  ae_list<T> * new_list = new ae_list<T>();
  
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
    ae_list_node<T> * node = _first;
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

template <typename T>
ae_list<T> * ae_list<T>::extract_ending_sublist( ae_list_node<T> * first_node )
{
  ae_list<T> * new_list = new ae_list<T>();
  
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
    ae_list_node<T> * node = first_node;
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

// Erase the list (nodes are deleted). Objects are deleted if delete_obj is true
template <typename T>
void ae_list<T>::erase( bool delete_obj )
{
  //  printf("nb elts in list %d\n", _nb_elts);

  ae_list_node<T>* node = _first;

  while ( node != NULL )
  {
    ae_list_node<T> * next = node->_next;

    if( delete_obj )
    {
      assert( node->_obj != NULL );
      delete node->_obj;
      node->_obj = NULL;
    }
    
    delete node;
    node = next;
  }

  _last = _first = NULL;
  _nb_elts = 0;
}

template <typename T>
void ae_list<T>::merge( ae_list<T> * append )
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

template <typename T>
inline void ae_list<T>::add_list( ae_list<T> * const to_add )
{
  if ( to_add == NULL ) return; // Nothing to do
  
  ae_list_node<T> * to_add_node = to_add->get_first();
  
  while ( to_add_node != NULL )
  {
    add( new ae_list_node<T>( *to_add_node ) );
    
    to_add_node = to_add_node->get_next();
  }
}

template <typename T>
bool ae_list<T>::is_empty( void ) const
{
  return ( _first == NULL );
}


template <typename T>
int32_t ae_list<T>::get_nb_elts( void ) const
{
  return _nb_elts;
}


/*!
  Returns the position of the object provided or -1 if not found
 */
template <typename T>
int32_t ae_list<T>::get_position( T* obj ) const
{
  int32_t i = 0;

  ae_list_node<T> * node = _first;
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
template <typename T>
ae_list_node<T> * ae_list<T>::get_node( int32_t pos ) const
{
  if ( pos >= _nb_elts ) return NULL;

  ae_list_node<T> * node = _first;
  for ( int32_t i = 0 ; i < pos ; i++ )
  {
    node = node->_next;
  }
  return node;
}


/*!
  Returns the object at position <pos> or NULL if pos is invalid
 */
template <typename T>
T* ae_list<T>::get_object( int32_t pos ) const
{
  if ( pos >= _nb_elts ) return NULL;

  ae_list_node<T> * node = _first;
  for ( int32_t i = 0 ; i < pos ; i++ )
  {
    node = node->_next;
  }
  return node->_obj;
}


/*!
  Returns a pointer to the first object in the list where needle was found
  Returns NULL if needle was not found
  comparator : function that returns true if <value> was found in <object>
 */
template <typename T>
ae_list_node<T>* ae_list<T>::bsearch( void* needle, int ( * comparator ) ( const void * value, const void * object ) ) const
{
  ae_list_node<T>* node = _first;

  while ( node != NULL )
  {
    if ( comparator( needle, (void*)node->get_obj() ) == 0 ) return node;

    node = node->_next;
  }

  return NULL;
}




//******************************************************************************
//                                  ae_list_node
//******************************************************************************
template<typename T>
inline ae_list_node<T>::ae_list_node( T obj, ae_list_node* prev, ae_list_node* next )
{
  _obj  = obj;
  _prev = prev;
  _next = next;
};

template<typename T>
ae_list_node<T>::ae_list_node( const ae_list_node<T> &model )
{
  _obj  = model._obj;
  _prev = NULL;
  _next = NULL;
};

template<typename T>
ae_list_node<T>::~ae_list_node( void )
{
}

template<typename T>
inline T& ae_list_node<T>::get_obj( void )
{
  return _obj;
}

template<typename T>
inline ae_list_node<T>* ae_list_node<T>::get_prev( void )
{
  return _prev;
}

template<typename T>
inline ae_list_node<T>* ae_list_node<T>::get_next( void )
{
  return _next;
}

template<typename T>
inline void ae_list_node<T>::set_obj( T obj )
{
  _obj = obj;
}

#endif // AE_LIST_H
