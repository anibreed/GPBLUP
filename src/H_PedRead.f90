module H_PedRead
  use M_IO_Ped, only: PEDHashNode, PEDHashTable, &
      pht_create, pht_insert_with_key, pht_search, pht_delete, pht_free, &
      pht_print_stats, pht_load_ped_file
  implicit none

  public :: PEDHashNode, PEDHashTable
  public :: pht_create, pht_insert_with_key, pht_search, pht_delete, pht_free
  public :: pht_print_stats, pht_load_ped_file
end module H_PedRead
