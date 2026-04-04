module M_readpar
  use M_readpar_relped, only: relped_params_t, read_relped_params
  use M_readpar_popqc, only: sentence, ParSections, NAMEKEY, read_parameters, read_to_sentences, &
      M_readpar_get_popqc_thresholds
  implicit none

    public :: relped_params_t, read_relped_params
  public :: read_parameters, read_to_sentences, sentence, ParSections, NAMEKEY, &
      M_readpar_get_popqc_thresholds
end module M_readpar
