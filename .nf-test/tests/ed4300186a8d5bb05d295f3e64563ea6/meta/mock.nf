import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { UMI_count } from '/kyukon/scratch/gent/vo/000/gvo00027/vsc44685/panoramaseq/nf-core-panoramaseq/modules/local/UMI_count/tests/../main.nf'

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
    
                input[0] = [
                    [ id:'S09', N_barcodes:34500, sample_size:10000, len_barcode:36, Nthresh:9, Ntriage:1000 ], // meta map
                    file('https://raw.githubusercontent.com/francops1722/test-datasets/panoramaseq/testdata/samtools/S09.assigned_sorted.bam', checkIfExists: true),
                    file('https://raw.githubusercontent.com/francops1722/test-datasets/panoramaseq/testdata/samtools/S09.assigned_sorted.bam.bai', checkIfExists: true)
                ]
                
    //----

    //run process
    UMI_count(*input)

    if (UMI_count.output){

        // consumes all named output channels and stores items in a json file
        for (def name in UMI_count.out.getNames()) {
            serializeChannel(name, UMI_count.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = UMI_count.out as Object[]
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
