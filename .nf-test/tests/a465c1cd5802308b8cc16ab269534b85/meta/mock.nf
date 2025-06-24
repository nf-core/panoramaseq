import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { UMITOOLS_EXTRACT } from '/kyukon/scratch/gent/vo/000/gvo00027/vsc44685/panoramaseq/nf-core-panoramaseq/modules/nf-core/umitools/extract/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    

    // process mapping
    def input = []
    
                input[0] = [ [ id:'test', single_end:false ], // meta map
                    [ file(params.modules_testdata_base_path + "genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true),
                    file(params.modules_testdata_base_path + "genomics/sarscov2/illumina/fastq/test_2.fastq.gz", checkIfExists: true) ]
                ]
                
    //----

    //run process
    UMITOOLS_EXTRACT(*input)

    if (UMITOOLS_EXTRACT.output){

        // consumes all named output channels and stores items in a json file
        for (def name in UMITOOLS_EXTRACT.out.getNames()) {
            serializeChannel(name, UMITOOLS_EXTRACT.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = UMITOOLS_EXTRACT.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
