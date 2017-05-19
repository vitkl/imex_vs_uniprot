<map version="1.0.1">
<!-- To view this file, download free mind mapping software FreeMind from http://freemind.sourceforge.net -->
<node CREATED="1480676247301" ID="ID_749205698" MODIFIED="1480676282850" TEXT="proteome_vs_interactome analysis">
<node CREATED="1480676294087" ID="ID_844643668" MODIFIED="1483546230970" POSITION="right">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>proteome_vs_IMEx(or IntAct) </b>
    </p>
    <p>
      <b>+ proteome_vs_IMEx_vs_Missing proteins (as annotated by Uniprot)</b>
    </p>
  </body>
</html></richcontent>
<node CREATED="1480679928831" ID="ID_1292174962" MODIFIED="1483546231178">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>intact vs uniprot analysis.R </b>
    </p>
    <p>
      script
    </p>
    <p>
      - calls other functions to download, clean,
    </p>
    <p>
      transform and analyze data
    </p>
    <p>
      - loops over: Species, Uniprot reviewed status,
    </p>
    <p>
      incl./excl. isoforms
    </p>
    <p>
      - requires to specify: date, missing/present protein evidence
    </p>
    <p>
      - plots venn diagrams and bar-plots
    </p>
    <p>
      - draft code: intact_vs_uniprot analysis - unused draft code.R
    </p>
  </body>
</html></richcontent>
<node CREATED="1480679931190" ID="ID_1532004159" MODIFIED="1481134725906">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>intact_vs_uniprot.R </b>
    </p>
    <p>
      (SPECIES_NAME , reviewed, isoforms, missing_proteins, date)
    </p>
    <p>
      can take only single value for each argument (vector length = 1)
    </p>
    <p>
      given vector length &gt; 1 produces weird results
    </p>
    <p>
      date is Sys.Date(), missing_proteins = TRUE by default
    </p>
  </body>
</html></richcontent>
<node CREATED="1480680510732" FOLDED="true" ID="ID_1009041512" MODIFIED="1481134735042">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>SPECIES_NAME_TO_ID </b>
    </p>
    <p>
      (SPECIES_NAME)
    </p>
    <p>
      finds SPECIES_ID and&#160;Proteome_ID given SPECIES_NAME
    </p>
  </body>
</html></richcontent>
<node CREATED="1480682145659" ID="ID_1450015653" MODIFIED="1480682162831">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      finds SPECIES_ID and&#160;Proteome_ID given SPECIES_NAME
    </p>
    <p>
      from:
    </p>
    <p>
      /uniprot/current_release/knowledgebase/reference_proteomes/README
    </p>
    <p>
      saves that file with the download date
    </p>
    <p>
      returns data.frame with SPECIES_ID and&#160;Proteome_ID
    </p>
    <p>
      usage:
    </p>
    <p>
      &#160;&#160;source(&quot;SPECIES_NAME_TO_ID.R&quot;)
    </p>
    <p>
      &#160;&#160;SPECIES_IDs = SPECIES_NAME_TO_ID(SPECIES_NAME)
    </p>
    <p>
      &#160;&#160;SPECIES_ID = SPECIES_IDs$SPECIES_ID;
    </p>
    <p>
      &#160;&#160;Proteome_ID = SPECIES_IDs$Proteome_ID;
    </p>
  </body>
</html></richcontent>
</node>
</node>
<node CREATED="1480680793711" ID="ID_1323367254" MODIFIED="1481134945303">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>query_PSICQUIC_for_interactions</b>
    </p>
    <p>
      (SPECIES_ID, SPECIES_NAME, databases, date)
    </p>
    <p>
      gets interactions for given species using PSICQUIC, saves&#160;and returns data.frame in MI-TAB 2.5 format
    </p>
  </body>
</html></richcontent>
<node CREATED="1480682279053" ID="ID_1349115463" MODIFIED="1480686393451">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      (SPECIES_ID = NA, SPECIES_NAME = &quot;Homo sapiens&quot;, databases = NA, date = Sys.Date())
    </p>
    <p>
      - if SPECIES_ID not provided calls SPECIES_NAME_TO_ID
    </p>
    <p>
      - queries &quot;databases&quot; for interactions for a SPECIES_ID using PSICQUIC&#160;package
    </p>
    <p>
      - rbinds query results from different databases
    </p>
    <p>
      - saves total query result into binary file (/Data) called:<br />&lt;databaseName_&quot;first database in list&quot;...&quot;last database in list&quot;_speciesID_&quot;SPECIES_ID&quot;_&quot;SPECIES_NAME&quot;_&quot;date&quot;&gt;
    </p>
    <p>
      - writes log (/Data/logs) file showing which databases yielded no interactions
    </p>
    <p>
      - prints the number of interactions/database to the prompt
    </p>
    <p>
      
    </p>
    <p>
      - before using PSICQUIC - checks if there is a saved binary file
    </p>
    <p>
      with interactions for given species, date, databases
    </p>
    <p>
      - if not - queries PSICQUIC
    </p>
    <p>
      
    </p>
    <p>
      - if yes - returns interactions in data.frame in format MI-TAB-2.5
    </p>
    <p>
      
    </p>
    <p>
      -has default database list (all IMEx) and default species (human), date = Sys.Date()
    </p>
    <p>
      
    </p>
    <p>
      usage:
    </p>
    <p>
      &#160;&#160;&#160;&#160;databases &lt;- c(&quot;mentha&quot;)
    </p>
    <p>
      &#160;&#160;&#160;&#160;source(&quot;query_PSICQUIC_for_interactions.R&quot;)
    </p>
    <p>
      &#160;&#160;&#160;&#160;all_interactions = query_PSICQUIC_for_interactions(SPECIES_ID = SPECIES_ID,
    </p>
    <p>
      &#160;&#160;&#160;&#160;SPECIES_NAME = SPECIES_NAME, databases = databases, date)
    </p>
  </body>
</html></richcontent>
</node>
</node>
<node CREATED="1480680673195" ID="ID_1494256375" MODIFIED="1480686140655">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>download_reference_proteome_q </b>
    </p>
    <p>
      (SPECIES_ID, Proteome_ID, date)
    </p>
    <p>
      downloads reference proteome(proteome:&quot;xxx&quot;) from Uniprot using url query
    </p>
  </body>
</html></richcontent>
<node CREATED="1480685211962" ID="ID_241886576" MODIFIED="1480686346492">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      (SPECIES_ID = &quot;9606&quot;, Proteome_ID = &quot;UP000005640&quot;, date = Sys.Date())
    </p>
    <p>
      - before sending url query - checks if there is a saved tab-delimited file with proteins for given proteome, date
    </p>
    <p>
      - if not - downloads reference proteome using url query
    </p>
    <p>
      (Uniprot search result table in tab format)
    </p>
    <p>
      - reads saved file into data.frame
    </p>
    <p>
      - prints number of proteins to the prompt
    </p>
    <p>
      
    </p>
    <p>
      - url query specifies proteome:&quot;Proteome_ID&quot;
    </p>
    <p>
      - url query specifies columns=id,reviewed,organism,proteome,organism-id,
    </p>
    <p>
      protein%20names,comment(ALTERNATIVE%20PRODUCTS),existence
    </p>
    <p>
      
    </p>
    <p>
      - arguments have default values for human, date = Sys.Date()
    </p>
    <p>
      
    </p>
    <p>
      usage:
    </p>
    <p>
      &#160;&#160;&#160;&#160;&#160;source(&quot;download_reference_proteome_q.R&quot;)
    </p>
    <p>
      &#160;&#160;&#160;&#160;&#160;reference_proteome_query = download_reference_proteome_q(SPECIES_ID, Proteome_ID, date)
    </p>
  </body>
</html></richcontent>
</node>
</node>
<node CREATED="1480685981322" ID="ID_1138326105" MODIFIED="1480686098166">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>download_whole_proteome </b>
    </p>
    <p>
      (SPECIES_ID, date)
    </p>
    <p>
      downloads whole proteome(species:&quot;xxx&quot;) from Uniprot using url query
    </p>
  </body>
</html></richcontent>
<node CREATED="1480686136673" ID="ID_635690856" MODIFIED="1480686584370">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      (SPECIES_ID = &quot;9606&quot;, date = Sys.Date())
    </p>
    <p>
      - before sending url query - checks if there is a saved tab-delimited file with proteins for given SPECIES_ID, date
    </p>
    <p>
      - if not - downloads reference proteome using url query
    </p>
    <p>
      (Uniprot search result table in tab format)
    </p>
    <p>
      - reads saved file into data.frame
    </p>
    <p>
      - prints number of proteins to the prompt
    </p>
    <p>
      
    </p>
    <p>
      - url query specifies species:&quot;SPECIES_ID&quot;
    </p>
    <p>
      - url query specifies columns
    </p>
    <p>
      =id,reviewed,organism,proteome,organism-id,
    </p>
    <p>
      protein%20names,comment(ALTERNATIVE%20PRODUCTS),existence
    </p>
    <p>
      
    </p>
    <p>
      - SPECIES_ID argument has default value for human, date&#160;= Sys.Date()
    </p>
    <p>
      
    </p>
    <p>
      usage:
    </p>
    <p>
      &#160;&#160;&#160;&#160;&#160;&#160;&#160;source(&quot;download_whole_proteome.R&quot;)
    </p>
    <p>
      &#160;&#160;&#160;&#160;&#160;&#160;&#160;whole_proteome_query = download_whole_proteome(SPECIES_ID, date)
    </p>
  </body>
</html></richcontent>
</node>
</node>
<node CREATED="1480686629290" ID="ID_1036966975" MODIFIED="1481134922675" TEXT="cleaning and filtering">
<node CREATED="1480686714838" ID="ID_454987379" MODIFIED="1480687622381">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      select all UniprotKB/SwissProt/trEMBL
    </p>
    <p>
      
    </p>
    <p>
      if reviewed = 2(SwissProt) or reviewed = 3(trEMBL)
    </p>
    <p>
      =&gt; proteome is filtered by SwissProt/trEMBL
    </p>
    <p>
      filtering will fail if SwissProt or trEMBL entries don't exist
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480686883617" ID="ID_1536711596" MODIFIED="1480687487170">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      isoforms - Uniprot - add or remove
    </p>
  </body>
</html></richcontent>
<node CREATED="1480686918968" FOLDED="true" ID="ID_911201573" MODIFIED="1480948545146">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      extract and add isoforms to proteome
    </p>
    <p>
      remove isoform-1
    </p>
  </body>
</html></richcontent>
<node CREATED="1480687221025" ID="ID_123140225" MODIFIED="1480688057499">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>isoform_id_extractor_Uniprot </b>
    </p>
    <p>
      (proteome)
    </p>
    <p>
      takes Uniprot query result table (data.frame)
    </p>
    <p>
      extracts isoforms from Uniprot column
    </p>
    <p>
      &quot;comment(ALTERNATIVE%20PRODUCTS)&quot;
    </p>
    <p>
      uses strapplyc
    </p>
    <p>
      and UniprotAC regular expression: &quot;IsoId=([[:alnum:]]+-[[:digit:]]+);&quot;
    </p>
    <p>
      
    </p>
    <p>
      usage:
    </p>
    <p>
      &#160;&#160;&#160;&#160;source(&quot;isoform_id_extractor_Uniprot.R&quot;)
    </p>
    <p>
      &#160;&#160;&#160;&#160;reference_proteome_isoforms = isoform_id_extractor_Uniprot(reference_proteome_query)
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480687267307" ID="ID_1509494625" MODIFIED="1480687414845">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      combine generic UniprotAC and isoforms UniprotAC
    </p>
    <p>
      in the character vector (other information from Uniprot search result is lost as this stage)
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480687511271" ID="ID_1596004615" MODIFIED="1480688726705">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>isoform_id_1_remover</b>
    </p>
    <p>
      (isoform_ids_1)
    </p>
    <p>
      replaces/removes isoform id-1 with generic: Q65483-1 =&gt; Q65483
    </p>
    <p>
      &#160;gsub(&quot;-1$&quot;,&quot;&quot;,isoform_ids_1)
    </p>
    <p>
      
    </p>
    <p>
      usage:
    </p>
    <p>
      &#160;&#160;&#160;&#160;source(&quot;isoform_id_1_remover.R&quot;)
    </p>
    <p>
      &#160;&#160;&#160;&#160;all_interactors_SPECIES_ID_only$interactor_IDs =
    </p>
    <p>
      &#160;&#160;&#160;&#160;isoform_id_1_remover(all_interactors_SPECIES_ID_only$interactor_IDs)
    </p>
    <p>
      &#160;&#160;&#160;&#160;reference_proteome = isoform_id_1_remover(reference_proteome)
    </p>
    <p>
      &#160;&#160;&#160;&#160;whole_proteome = isoform_id_1_remover(whole_proteome)
    </p>
  </body>
</html></richcontent>
</node>
</node>
</node>
<node CREATED="1480687624921" ID="ID_1825843302" MODIFIED="1480688309976">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>interactions_to_interactors </b>
    </p>
    <p>
      (all_interactions)
    </p>
    <p>
      transform the table of interactions into the table of interactors
    </p>
    <p>
      cleans information from MI-TAB labels
    </p>
    <p>
      saves the following information:<br />interactor_IDs - interactor identifier
    </p>
    <p>
      interactor_IDs_databases - type of interactor identifier
    </p>
    <p>
      interactor_databases - database from where interactor is coming
    </p>
    <p>
      interactor_SPECIES_ID - SPECIES_ID of interactor from interactor_databases
    </p>
    <p>
      
    </p>
    <p>
      usage:
    </p>
    <p>
      &#160;&#160;source(&quot;interactions_to_interactors.R&quot;)
    </p>
    <p>
      &#160;&#160;all_interactors = interactions_to_interactors(all_interactions)
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480688314575" ID="ID_333385792" MODIFIED="1480688534695">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>uniprotkb_and_SPECIES_ID_interactor_selector</b>
    </p>
    <p>
      (all_interactors, SPECIES_ID)
    </p>
    <p>
      filters interactor list from identifiers which are not UniprotKB and do not belong to given species
    </p>
    <p>
      uses dplyr::filter
    </p>
    <p>
      prints summary of how many interactors have which identifiers to the prompt
    </p>
    <p>
      
    </p>
    <p>
      usage:
    </p>
    <p>
      &#160;&#160;source(&quot;uniprotkb_and_SPECIES_ID_interactor_selector.R&quot;)
    </p>
    <p>
      &#160;&#160;all_interactors_SPECIES_ID_only = uniprotkb_and_SPECIES_ID_interactor_selector(all_interactors, SPECIES_ID)
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480688551220" ID="ID_1780160808" MODIFIED="1480688753437">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      if isoforms are included in the analysis:
    </p>
    <p>
      isoform-1 IDs are replaced with generic (Uniprot) interactor identifiers
    </p>
    <p>
      if isoforms are NOT included in the analysis
    </p>
    <p>
      all isoformIDs are replaced with generic
    </p>
  </body>
</html></richcontent>
<node CREATED="1480688730988" ID="ID_1943744081" MODIFIED="1480688749941">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>isoform_id_1_remover&#160;</b>
    </p>
    <p>
      (isoform_ids_1)
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480688754769" ID="ID_1615882632" MODIFIED="1480688909367">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>isoform_id_all_remover </b>
    </p>
    <p>
      (isoform_ids)
    </p>
    <p>
      replaces/removes isoform all id with generic: Q65483-1 =&gt; Q65483
    </p>
    <p>
      gsub(&quot;-[[:digit:]]+$&quot;,&quot;&quot;,isoform_ids)
    </p>
    <p>
      
    </p>
    <p>
      usage:
    </p>
    <p>
      &#160;&#160;if(isoforms == FALSE){
    </p>
    <p>
      &#160;&#160;&#160;&#160;source(&quot;isoform_id_all_remover.R&quot;)
    </p>
    <p>
      &#160;&#160;&#160;&#160;all_interactors_SPECIES_ID_only$interactor_IDs = isoform_id_all_remover(all_interactors_SPECIES_ID_only$interactor_IDs)
    </p>
    <p>
      &#160;&#160;}
    </p>
  </body>
</html></richcontent>
</node>
</node>
</node>
<node CREATED="1480689282746" ID="ID_1678171279" MODIFIED="1480689906008">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      Transforming data into logic table
    </p>
    <p>
      (presence/absence of protein X in database (or list) Y = 1/0)
    </p>
  </body>
</html></richcontent>
<node CREATED="1480689406360" ID="ID_796049181" MODIFIED="1480689897646">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>proteome_vs_interactome_logic_table_generator</b>
    </p>
    <p>
      
    </p>
    <p>
      uses tidyr::spread to transform
    </p>
    <p>
      P86453&#160;&#160;&#160;&#160;IntAct&#160;&#160;&#160;&#160;1
    </p>
    <p>
      P42789&#160;&#160;&#160;&#160;MINT&#160;&#160;&#160;&#160;&#160;1
    </p>
    <p>
      P83230&#160;&#160;&#160;&#160;IntAct&#160;&#160;&#160;&#160;1
    </p>
    <p>
      
    </p>
    <p>
      into
    </p>
    <p>
      &#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;IntAct&#160;&#160;&#160;&#160;MINT
    </p>
    <p>
      P86453&#160;&#160;&#160;&#160;1&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;0
    </p>
    <p>
      P42789&#160;&#160;&#160;&#160;0&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;1
    </p>
    <p>
      P83230&#160;&#160;&#160;&#160;1&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;&#160;0
    </p>
    <p>
      (tidyr::spread allows to specify to fill empty intersects with 0-s)
    </p>
    <p>
      
    </p>
    <p>
      then merges interactome table to proteome table using merge
    </p>
    <p>
      (merge fills empty intersects with NA-s)
    </p>
    <p>
      then changes NA generated by merge into 0-s
    </p>
    <p>
      returns data.frame
    </p>
    <p>
      
    </p>
    <p>
      usage:
    </p>
    <p>
      &#160;&#160;&#160;&#160;source(&quot;proteome_vs_interactome_logic_table_generator.R&quot;)
    </p>
    <p>
      &#160;&#160;&#160;&#160;proteome_vs_interactome_f = proteome_vs_interactome_logic_table_generator(reference_proteome_unique,
    </p>
    <p>
      &#160;&#160;&#160;&#160;whole_proteome_unique, all_interactors_SPECIES_ID_only, unique_interactors_SPECIES_ID_only)
    </p>
  </body>
</html></richcontent>
</node>
</node>
<node CREATED="1480688917193" ID="ID_193641797" MODIFIED="1480688926140" TEXT="missing protein evidence">
<node CREATED="1480689116189" ID="ID_1350604556" MODIFIED="1480689181726" TEXT="to create a missing protein list">
<node CREATED="1480688959452" ID="ID_796802181" MODIFIED="1481134936874" TEXT="if(missing_proteins == FALSE)">
<node CREATED="1480688963411" ID="ID_647706995" MODIFIED="1480689114649" TEXT="proteome from query is filtered by the presence of Evidence at protein level"/>
</node>
<node CREATED="1480689094052" ID="ID_812787235" MODIFIED="1480689103982" TEXT="if(missing_proteins == TRUE)">
<node CREATED="1480689118933" ID="ID_1694052738" MODIFIED="1480689130117">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      proteome from query is filtered by the absence of Evidence at protein level
    </p>
  </body>
</html></richcontent>
</node>
</node>
</node>
<node CREATED="1480689240320" ID="ID_1601053707" MODIFIED="1480689961818">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      merge missing protein list with the total logic table
    </p>
    <p>
      using merge function
    </p>
    <p>
      then change NA generated by merge into 0-s
    </p>
  </body>
</html></richcontent>
</node>
</node>
<node CREATED="1480690051828" ID="ID_1412062482" MODIFIED="1480690210672">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      save logic table as a tab delimited file
    </p>
    <p>
      
    </p>
    <p>
      will be unique for each value of SPECIES_NAME, reviewed, isoforms, date
    </p>
    <p>
      
    </p>
    <p>
      file is called
    </p>
    <p>
      proteome_vs_interactome_f_&quot;SPECIES_ID&quot;_reviewed_&quot;1/2/3&quot;_isoforms_&quot;TRUE/FALSE&quot;_&quot;date&quot;.txt
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480690232350" ID="ID_978996028" MODIFIED="1480690378039" TEXT="print and save summary of how many interactors have non-Uniprot identifiers and how many interactors come not from given species">
<node CREATED="1480690314800" ID="ID_1748848970" MODIFIED="1480690350602">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>uniprotkb_and_SPECIES_ID_interactor_summary</b>
    </p>
    <p>
      (all_interactors, SPECIES_ID, SPECIES_NAME)
    </p>
    <p>
      
    </p>
    <p>
      &#160;&#160;&#160;&#160;source(&quot;uniprotkb_and_SPECIES_ID_interactor_summary.R&quot;)
    </p>
    <p>
      &#160;&#160;&#160;&#160;uniprotkb_and_SPECIES_ID_interactor_summary = uniprotkb_and_SPECIES_ID_interactor_summary(all_interactors, SPECIES_ID, SPECIES_NAME)
    </p>
  </body>
</html></richcontent>
</node>
</node>
</node>
<node CREATED="1480680450316" FOLDED="true" ID="ID_1554778389" MODIFIED="1480932771935">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>intact_vs_uniprot_overlap.R</b>
    </p>
    <p>
      (SPECIES_NAME, SPECIES_ID, reviewed, isoforms)
    </p>
    <p>
      Does overlap comparisons and saves the result
    </p>
    <p>
      
    </p>
    <p>
      -reads logic table generated by intact_vs_uniprot.R
    </p>
  </body>
</html></richcontent>
<node CREATED="1480690537336" ID="ID_706331853" MODIFIED="1480690665356">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      (SPECIES_NAME, SPECIES_ID, reviewed, isoforms)
    </p>
    <p>
      
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480690713904" FOLDED="true" ID="ID_1320082877" MODIFIED="1480932766986">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      proteome_vs_IMEx(or IntAct)
    </p>
  </body>
</html></richcontent>
<node CREATED="1480690632758" ID="ID_1769177745" MODIFIED="1480692314752">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>proteome_vs_interactome_summary </b>
    </p>
    <p>
      (proteome_vs_interactome_f, database = &quot;IMEx&quot;, SPECIES_NAME, reviewed, isoforms)
    </p>
    <p>
      
    </p>
    <p>
      does pairwise overlap comparisons
    </p>
    <p>
      old function - works only to compare whole_proteome_Uniprot and reference_proteome_Uniprot with database provided as an argument = &quot;IMEx&quot;
    </p>
  </body>
</html></richcontent>
</node>
</node>
<node CREATED="1480690776823" FOLDED="true" ID="ID_535115111" MODIFIED="1480932768911" TEXT="proteome_vs_IMEx_vs_Missing proteins (as annotated by Uniprot) ">
<node CREATED="1480690642421" ID="ID_1609416996" MODIFIED="1480690656333">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>A_vs_B_vs_C_overlap</b>
    </p>
  </body>
</html></richcontent>
</node>
</node>
</node>
<node CREATED="1480680410432" ID="ID_1723537875" MODIFIED="1480680410432" TEXT=""/>
</node>
</node>
<node CREATED="1480676381469" ID="ID_1304043330" MODIFIED="1480934525964" POSITION="right">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>proteome_vs_IMEx_vs_Biogrid </b>
    </p>
    <p>
      requires &quot;logic table&quot; from proteome_vs_IMEx
    </p>
    <p>
      &quot;proteome_vs_interactome_f_&quot;
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480934836767" ID="ID_629971211" MODIFIED="1480934912947" POSITION="right">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>proteome_vs_IMEx_BIND </b>or <b>proteome_vs_IMEx_HPRD </b>
    </p>
    <p>
      pending
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480950084920" ID="ID_265114987" MODIFIED="1480950284920" POSITION="right">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>missing_evidence_vs_IMEx_statistics</b>
    </p>
    <p>
      are proteins without protein evidence likely not to appear in IMEx?
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1480955468453" ID="ID_1206991431" MODIFIED="1481191476119" POSITION="right">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>swissprot_vs_imex_protein_properties </b>
    </p>
    <p>
      mass
    </p>
    <p>
      number of natural variants
    </p>
    <p>
      number of isoforms
    </p>
    <p>
      annotation score
    </p>
    <p>
      requires logic table from <b>proteome_vs_IMEx </b>(E.coli)<b>, </b>
    </p>
    <p>
      <b>proteome_vs_IMEx_vs_Biogrid </b>(other species)
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1481191379202" ID="ID_698080164" MODIFIED="1481191428624" POSITION="right">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>swissprot_vs_imex_interaction_properties </b>
    </p>
    <p>
      requires logic table from <b>swissprot_vs_imex_protein_properties</b>
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1481716218919" ID="ID_1063533877" MODIFIED="1481716247568" POSITION="right">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>swissprot_vs_imex_interaction_properties_bait_prey</b>
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1481716249404" ID="ID_760523772" MODIFIED="1481716269040" POSITION="right">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      <b>swissprot_vs_imex_interaction_properties_GO</b>
    </p>
  </body>
</html></richcontent>
</node>
<node CREATED="1481716270508" ID="ID_1587874182" MODIFIED="1481716283764" POSITION="right">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      swissprot_vs_imex_publications
    </p>
  </body>
</html></richcontent>
</node>
</node>
</map>
