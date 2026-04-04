module H_MapRead
  use M_IO_Map, only: MAPHashNode, MAPHashTable, &
      mht_create, mht_insert, mht_search, mht_free, mht_load_map_file
  implicit none

  public :: MAPHashNode, MAPHashTable
  public :: mht_create, mht_insert, mht_search, mht_free, mht_load_map_file
end module H_MapRead
